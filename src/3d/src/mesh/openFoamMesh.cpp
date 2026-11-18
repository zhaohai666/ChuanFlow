#include "cf/fvm/OpenFoamMesh.h"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "cf/common/Error.h"
#include "cf/common/Log.h"

namespace cf
{
namespace
{

// ---------------------------------------------------------------------------
// Small lexer for OpenFOAM's ASCII list syntax. Words, parentheses, braces and
// semicolons are tokens; C++-style comments and all whitespace are dropped.
// ---------------------------------------------------------------------------

struct Token
{
    enum Kind { Word, LParen, RParen, LBrace, RBrace, Semi };

    Kind kind = Word;
    std::string text;
};

bool isSpace(char c)
{
    return std::isspace(static_cast<unsigned char>(c)) != 0;
}

bool isDelimiter(char c)
{
    if (isSpace(c)) return true;
    return c == '(' || c == ')' || c == '{' || c == '}' || c == ';';
}

bool commentStart(const std::string& s, std::size_t i)
{
    return i + 1 < s.size() && s[i] == '/' && (s[i + 1] == '/' || s[i + 1] == '*');
}

std::vector<Token> lex(const std::string& text, const std::string& where)
{
    std::vector<Token> out;
    std::size_t i = 0;
    while (i < text.size())
    {
        if (commentStart(text, i))
        {
            if (text[i + 1] == '/')
            {
                while (i < text.size() && text[i] != '\n') ++i;
            }
            else
            {
                const std::size_t end = text.find("*/", i + 2);
                if (end == std::string::npos)
                    cfFatal("openfoam: unterminated /* comment in " + where);
                i = end + 2;
            }
            continue;
        }

        const char c = text[i];
        if (isSpace(c)) { ++i; continue; }

        const Token::Kind punct =
            c == '(' ? Token::LParen :
            c == ')' ? Token::RParen :
            c == '{' ? Token::LBrace :
            c == '}' ? Token::RBrace :
            c == ';' ? Token::Semi : Token::Word;
        if (punct != Token::Word)
        {
            out.push_back({punct, std::string(1, c)});
            ++i;
            continue;
        }

        std::size_t j = i;
        while (j < text.size() && !isDelimiter(text[j]) && !commentStart(text, j)) ++j;
        if (j == i) ++j; // defensive: never loop forever on a stray byte
        out.push_back({Token::Word, text.substr(i, j - i)});
        i = j;
    }
    return out;
}

std::string lower(std::string s)
{
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

std::string unquote(const std::string& s)
{
    if (s.size() >= 2 && (s.front() == '"' || s.front() == '\'')
        && s.back() == s.front())
    {
        return s.substr(1, s.size() - 2);
    }
    return s;
}

bool isIntegerWord(const std::string& s)
{
    if (s.empty()) return false;
    std::size_t i = (s[0] == '+' || s[0] == '-') ? 1 : 0;
    if (i == s.size()) return false;
    for (; i < s.size(); ++i)
    {
        if (!std::isdigit(static_cast<unsigned char>(s[i]))) return false;
    }
    return true;
}

long parseInt(const std::string& tok, const std::string& where)
{
    try
    {
        std::size_t used = 0;
        const long v = std::stol(tok, &used);
        if (used != tok.size()) cfFatal("openfoam: bad integer '" + tok + "' in " + where);
        return v;
    }
    catch (const std::exception&)
    {
        cfFatal("openfoam: bad integer '" + tok + "' in " + where);
    }
    return 0;
}

scalar parseScalar(const std::string& tok, const std::string& where)
{
    try
    {
        std::size_t used = 0;
        const double v = std::stod(tok, &used);
        if (used != tok.size()) cfFatal("openfoam: bad number '" + tok + "' in " + where);
        return static_cast<scalar>(v);
    }
    catch (const std::exception&)
    {
        cfFatal("openfoam: bad number '" + tok + "' in " + where);
    }
    return 0.0;
}

std::string readFile(const std::string& path)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) cfFatal("openfoam: cannot open " + path);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

std::size_t matching(const std::vector<Token>& t, std::size_t open,
                     Token::Kind opening, Token::Kind closing,
                     const std::string& where)
{
    int depth = 0;
    for (std::size_t i = open; i < t.size(); ++i)
    {
        if (t[i].kind == opening) ++depth;
        else if (t[i].kind == closing)
        {
            if (--depth == 0) return i;
        }
    }
    cfFatal("openfoam: unbalanced delimiters in " + where);
    return 0;
}

/// Locate the `N (` that opens the single top-level list of a file. Tokens
/// inside the FoamFile header braces are skipped, so `version 2.0;` etc. can
/// never be mistaken for the count.
struct ListRegion
{
    long count = 0;
    std::size_t begin = 0; // first token inside the outer parentheses
    std::size_t end = 0;   // one past the last token inside
};

ListRegion findList(const std::vector<Token>& t, const std::string& where)
{
    int depth = 0;
    for (std::size_t i = 0; i + 1 < t.size(); ++i)
    {
        if (t[i].kind == Token::LBrace) ++depth;
        else if (t[i].kind == Token::RBrace) --depth;
        else if (depth == 0 && t[i].kind == Token::Word && isIntegerWord(t[i].text)
                 && t[i + 1].kind == Token::LParen)
        {
            const long n = parseInt(t[i].text, where);
            if (n < 0) cfFatal("openfoam: negative list size in " + where);

            const std::size_t close = matching(t, i + 1, Token::LParen, Token::RParen, where);
            ListRegion r;
            r.count = n;
            r.begin = i + 2;
            r.end = close;
            return r;
        }
    }
    cfFatal("openfoam: no 'N ( ... )' list in " + where);
    return {};
}

class Cursor
{
public:
    Cursor(const std::vector<Token>& t, std::size_t pos, std::size_t end, std::string where)
        : t_(t), pos_(pos), end_(end), where_(std::move(where))
    {
    }

    bool atEnd() const { return pos_ >= end_; }

    const Token& peek() const
    {
        if (atEnd()) cfFatal("openfoam: truncated list in " + where_);
        return t_[pos_];
    }

    Token take()
    {
        const Token tok = peek();
        ++pos_;
        return tok;
    }

    void expect(Token::Kind kind, const std::string& what)
    {
        if (peek().kind != kind)
            cfFatal("openfoam: expected " + what + " in " + where_);
        ++pos_;
    }

    long takeInt(const std::string& what)
    {
        const Token tok = take();
        if (tok.kind != Token::Word || !isIntegerWord(tok.text))
            cfFatal("openfoam: expected " + what + " in " + where_);
        return parseInt(tok.text, where_);
    }

    scalar takeScalar(const std::string& what)
    {
        const Token tok = take();
        if (tok.kind != Token::Word)
            cfFatal("openfoam: expected " + what + " in " + where_);
        return parseScalar(tok.text, where_);
    }

    std::size_t pos() const { return pos_; }
    void seek(std::size_t p) { pos_ = p; }

private:
    const std::vector<Token>& t_;
    std::size_t pos_;
    std::size_t end_;
    std::string where_;
};

// ---------------------------------------------------------------------------
// The five files
// ---------------------------------------------------------------------------

struct ParsedList
{
    std::vector<Token> tokens;
    ListRegion region;
    std::string where;
};

ParsedList loadList(const std::string& path)
{
    ParsedList out;
    out.where = path;
    out.tokens = lex(readFile(path), path);
    out.region = findList(out.tokens, path);
    return out;
}

std::vector<Vector> parsePoints(const ParsedList& data)
{
    Cursor c(data.tokens, data.region.begin, data.region.end, data.where);
    std::vector<Vector> pts;
    pts.reserve(static_cast<std::size_t>(data.region.count));
    for (long i = 0; i < data.region.count; ++i)
    {
        c.expect(Token::LParen, "point '('");
        const scalar x = c.takeScalar("point x");
        const scalar y = c.takeScalar("point y");
        const scalar z = c.takeScalar("point z");
        c.expect(Token::RParen, "point ')'");
        pts.push_back(Vector{x, y, z});
    }
    if (!c.atEnd()) cfFatal("openfoam: extra entries in points of " + data.where);
    return pts;
}

std::vector<std::vector<label>> parseFaces(const ParsedList& data)
{
    Cursor c(data.tokens, data.region.begin, data.region.end, data.where);
    std::vector<std::vector<label>> faces;
    faces.reserve(static_cast<std::size_t>(data.region.count));
    for (long i = 0; i < data.region.count; ++i)
    {
        const long n = c.takeInt("face vertex count");
        if (n < 3) cfFatal("openfoam: face " + std::to_string(i) + " has fewer than 3 vertices");
        c.expect(Token::LParen, "face '('");
        std::vector<label> face;
        face.reserve(static_cast<std::size_t>(n));
        for (long k = 0; k < n; ++k)
        {
            const long v = c.takeInt("face vertex index");
            if (v < 0) cfFatal("openfoam: negative point index in face " + std::to_string(i));
            face.push_back(static_cast<label>(v));
        }
        c.expect(Token::RParen, "face ')'");
        faces.push_back(std::move(face));
    }
    if (!c.atEnd()) cfFatal("openfoam: extra entries in faces of " + data.where);
    return faces;
}

std::vector<label> parseLabels(const ParsedList& data, const std::string& what)
{
    Cursor c(data.tokens, data.region.begin, data.region.end, data.where);
    std::vector<label> values;
    values.reserve(static_cast<std::size_t>(data.region.count));
    for (long i = 0; i < data.region.count; ++i)
    {
        const long v = c.takeInt(what);
        if (v < 0) cfFatal("openfoam: negative " + what + " in " + data.where);
        values.push_back(static_cast<label>(v));
    }
    if (!c.atEnd()) cfFatal("openfoam: extra entries in " + what + " of " + data.where);
    return values;
}

struct PatchSpec
{
    std::string name;
    std::string type;
    long nFaces = 0;
    long startFace = 0;
};

std::vector<PatchSpec> parseBoundary(const ParsedList& data)
{
    Cursor c(data.tokens, data.region.begin, data.region.end, data.where);
    std::vector<PatchSpec> patches;
    patches.reserve(static_cast<std::size_t>(data.region.count));

    for (long i = 0; i < data.region.count; ++i)
    {
        const Token name = c.take();
        if (name.kind != Token::Word)
            cfFatal("openfoam: expected patch name in " + data.where);
        c.expect(Token::LBrace, "patch '{'");

        const std::size_t close = matching(data.tokens, c.pos() - 1,
                                           Token::LBrace, Token::RBrace, data.where);
        PatchSpec p;
        p.name = unquote(name.text);
        bool haveType = false, haveCount = false, haveStart = false;

        // Scan key/value pairs loosely: only the three keys we care about are
        // interpreted, so extra entries such as `inGroups 1(wall);` are fine.
        for (std::size_t j = c.pos(); j < close; ++j)
        {
            if (data.tokens[j].kind != Token::Word || j + 1 >= close) continue;
            const std::string key = lower(data.tokens[j].text);
            const Token& value = data.tokens[j + 1];
            if (key == "type" && value.kind == Token::Word)
            {
                p.type = unquote(value.text);
                haveType = true;
            }
            else if (key == "nfaces" && value.kind == Token::Word)
            {
                p.nFaces = parseInt(value.text, data.where);
                haveCount = true;
            }
            else if (key == "startface" && value.kind == Token::Word)
            {
                p.startFace = parseInt(value.text, data.where);
                haveStart = true;
            }
        }
        if (!haveType || !haveCount || !haveStart)
        {
            cfFatal("openfoam: patch '" + p.name + "' needs type, nFaces and startFace");
        }
        c.seek(close + 1);
        patches.push_back(std::move(p));
    }
    if (!c.atEnd()) cfFatal("openfoam: extra entries in boundary of " + data.where);
    return patches;
}

std::string mapPatchType(const std::string& raw)
{
    const std::string key = lower(unquote(raw));
    if (key == "wall" || key == "67" || key == "patch") return "wall";
    if (key == "empty" || key == "symmetryplane" || key == "symmetry") return "symmetry";
    return unquote(raw); // conservative: do not invent a boundary condition
}

polyMesh buildMesh(const std::vector<Vector>& points,
                   std::vector<std::vector<label>> faces,
                   std::vector<label> owner,
                   std::vector<label> neighbour,
                   const std::vector<PatchSpec>& patches,
                   const std::map<std::string, std::string>& patchTypes,
                   const std::string& where)
{
    const label nPoints = static_cast<label>(points.size());
    const label nFaces = static_cast<label>(faces.size());

    if (static_cast<label>(owner.size()) != nFaces)
        cfFatal("openfoam: owner has " + std::to_string(owner.size())
                + " entries but faces has " + std::to_string(nFaces));

    for (label f = 0; f < nFaces; ++f)
    {
        for (const label v : faces[static_cast<std::size_t>(f)])
        {
            if (v < 0 || v >= nPoints)
                cfFatal("openfoam: face " + std::to_string(f) + " refers to point "
                        + std::to_string(v) + ", outside 0.." + std::to_string(nPoints - 1));
        }
    }

    label nCells = 0;
    for (const label o : owner) nCells = std::max(nCells, o + 1);
    for (const label nb : neighbour) nCells = std::max(nCells, nb + 1);
    if (nCells <= 0) cfFatal("openfoam: mesh " + where + " has no cells");

    std::vector<char> ownedCell(static_cast<std::size_t>(nCells), 0);
    for (const label o : owner)
    {
        if (o >= nCells)
            cfFatal("openfoam: owner index " + std::to_string(o) + " exceeds cell count");
        ownedCell[static_cast<std::size_t>(o)] = 1;
    }
    for (label c = 0; c < nCells; ++c)
    {
        if (!ownedCell[static_cast<std::size_t>(c)])
        {
            cfFatal("openfoam: cell " + std::to_string(c) + " owns no face"
                    " (owner index out of range?)");
        }
    }
    for (const label nb : neighbour)
    {
        if (nb < 0 || nb >= nCells)
            cfFatal("openfoam: neighbour index " + std::to_string(nb) + " exceeds cell count");
    }
    for (std::size_t f = 0; f < neighbour.size(); ++f)
    {
        if (owner[f] == neighbour[f])
            cfFatal("openfoam: face " + std::to_string(f) + " has owner == neighbour");
    }

    const label nInternal = static_cast<label>(neighbour.size());
    label nBoundary = 0;
    for (const PatchSpec& p : patches)
    {
        if (p.nFaces < 0 || p.startFace < 0 || p.startFace + p.nFaces > nFaces)
            cfFatal("openfoam: patch '" + p.name + "' range is outside the face list");
        nBoundary += p.nFaces;
    }
    if (nInternal + nBoundary != nFaces)
    {
        cfFatal("openfoam: neighbour (" + std::to_string(nInternal) + ") + boundary ("
                + std::to_string(nBoundary) + ") != faces (" + std::to_string(nFaces) + ")");
    }

    std::vector<char> covered(static_cast<std::size_t>(nFaces), 0);
    std::vector<polyMesh::Patch> meshPatches;
    meshPatches.reserve(patches.size());
    for (const PatchSpec& p : patches)
    {
        if (p.startFace < nInternal)
            cfFatal("openfoam: patch '" + p.name + "' covers internal face "
                    + std::to_string(p.startFace));
        for (label f = p.startFace; f < p.startFace + p.nFaces; ++f)
        {
            if (covered[static_cast<std::size_t>(f)])
                cfFatal("openfoam: boundary face " + std::to_string(f) + " is in two patches");
            covered[static_cast<std::size_t>(f)] = 1;
        }

        polyMesh::Patch mp;
        mp.name = p.name;
        const auto it = patchTypes.find(p.name);
        mp.type = it != patchTypes.end() ? it->second : mapPatchType(p.type);
        mp.start = p.startFace;
        mp.size = p.nFaces;
        meshPatches.push_back(std::move(mp));
    }
    for (label f = nInternal; f < nFaces; ++f)
    {
        if (!covered[static_cast<std::size_t>(f)])
            cfFatal("openfoam: boundary face " + std::to_string(f) + " is in no patch");
    }

    Field<Vector> pointsField(points.begin(), points.end());
    Field<label> ownerField(owner.begin(), owner.end());
    Field<label> neighbourField(neighbour.begin(), neighbour.end());

    polyMesh mesh(std::move(pointsField), std::move(faces), std::move(ownerField),
                  std::move(neighbourField), std::move(meshPatches));

    scalar volume = 0;
    for (const scalar v : mesh.cellVolumes()) volume += v;
    log::info("openfoam: ", where, " : cells ", mesh.nCells(), ", faces ", mesh.nFaces(),
              " (internal ", mesh.nInternalFaces(), ", boundary ", mesh.nBoundaryFaces(),
              "), points ", mesh.nPoints(), ", patches ", mesh.nPatches(),
              ", volume ", volume);
    return mesh;
}

} // namespace

polyMesh OpenFoamMesh::fromFiles(const std::string& pointsFile,
                                 const std::string& facesFile,
                                 const std::string& ownerFile,
                                 const std::string& neighbourFile,
                                 const std::string& boundaryFile,
                                 const std::map<std::string, std::string>& patchTypes)
{
    // Reject binary meshes early with an actionable message.
    for (const std::string& path : {pointsFile, facesFile, ownerFile, neighbourFile, boundaryFile})
    {
        const std::string text = readFile(path);
        const std::size_t fmt = lower(text).find("format");
        if (fmt != std::string::npos)
        {
            const std::size_t semi = text.find(';', fmt);
            if (semi != std::string::npos && lower(text.substr(fmt, semi - fmt)).find("binary")
                != std::string::npos)
            {
                cfFatal("openfoam: " + path + " is a binary mesh; write it with ASCII format");
            }
        }
    }

    const ParsedList points = loadList(pointsFile);
    const ParsedList faces = loadList(facesFile);
    const ParsedList owner = loadList(ownerFile);
    const ParsedList neighbour = loadList(neighbourFile);
    const ParsedList boundary = loadList(boundaryFile);

    return buildMesh(parsePoints(points), parseFaces(faces),
                     parseLabels(owner, "owner"), parseLabels(neighbour, "neighbour"),
                     parseBoundary(boundary), patchTypes, pointsFile);
}

polyMesh OpenFoamMesh::fromDirectory(const std::string& directory,
                                     const std::map<std::string, std::string>& patchTypes)
{
    namespace fs = std::filesystem;

    fs::path dir(directory);
    if (!fs::exists(dir / "points") && fs::exists(dir / "constant" / "polyMesh" / "points"))
    {
        dir = dir / "constant" / "polyMesh";
    }
    if (!fs::exists(dir / "points"))
    {
        cfFatal("openfoam: no 'points' file in '" + directory
                + "' (nor in its constant/polyMesh subdirectory)");
    }

    return fromFiles((dir / "points").string(), (dir / "faces").string(),
                     (dir / "owner").string(), (dir / "neighbour").string(),
                     (dir / "boundary").string(), patchTypes);
}

} // namespace cf
