#include "cf/common/Dictionary.h"

#include <cctype>
#include <cstdlib>
#include <fstream>
#include <sstream>

namespace cf
{
namespace
{

enum class TokType { word, number, lbrace, rbrace, semi, lparen, rparen, eof };

struct Token
{
    TokType type = TokType::eof;
    std::string text;
    double number = 0;
};

class Lexer
{
public:
    explicit Lexer(const std::string& src) : s_(src) {}

    Token next()
    {
        skipSpace();
        if (i_ >= s_.size())
        {
            return {TokType::eof, "", 0};
        }
        const char c = s_[i_];
        auto punct = [&](TokType t) { ++i_; return Token{t, std::string(1, c), 0}; };
        if (c == '{') return punct(TokType::lbrace);
        if (c == '}') return punct(TokType::rbrace);
        if (c == ';') return punct(TokType::semi);
        if (c == '(') return punct(TokType::lparen);
        if (c == ')') return punct(TokType::rparen);

        if (c == '"')
        {
            ++i_;
            std::string out;
            while (i_ < s_.size() && s_[i_] != '"')
            {
                out.push_back(s_[i_++]);
            }
            if (i_ < s_.size()) ++i_;
            return {TokType::word, out, 0};
        }

        std::string tok;
        while (i_ < s_.size() && !std::isspace(static_cast<unsigned char>(s_[i_]))
               && s_[i_] != '{' && s_[i_] != '}' && s_[i_] != ';'
               && s_[i_] != '(' && s_[i_] != ')')
        {
            tok.push_back(s_[i_++]);
        }

        char* end = nullptr;
        const double v = std::strtod(tok.c_str(), &end);
        if (!tok.empty() && end == tok.c_str() + tok.size())
        {
            return {TokType::number, tok, v};
        }
        return {TokType::word, tok, 0};
    }

private:
    void skipSpace()
    {
        for (;;)
        {
            while (i_ < s_.size() && std::isspace(static_cast<unsigned char>(s_[i_]))) ++i_;
            if (i_ < s_.size() && s_[i_] == '#')
            {
                while (i_ < s_.size() && s_[i_] != '\n') ++i_;
                continue;
            }
            if (i_ + 1 < s_.size() && s_[i_] == '/' && s_[i_ + 1] == '/')
            {
                while (i_ < s_.size() && s_[i_] != '\n') ++i_;
                continue;
            }
            return;
        }
    }

    const std::string& s_;
    std::size_t i_ = 0;
};

class Parser
{
public:
    explicit Parser(const std::string& src) : lex_(src) {}

    Dictionary parseTop()
    {
        Dictionary dict;
        for (;;)
        {
            Token key = lex_.next();
            if (key.type == TokType::eof) return dict;
            if (key.type != TokType::word)
            {
                cfFatal("dictionary: expected key, got '" + key.text + "'");
            }
            parseEntry(dict, key.text);
        }
    }

private:
    void parseEntry(Dictionary& dict, const std::string& key)
    {
        Dictionary::Node node;
        Token t = lex_.next();
        switch (t.type)
        {
            case TokType::number:
            {
                node.kind = Dictionary::Kind::scalar;
                node.scalar = t.number;
                node.word = t.text;
                expectSemi();
                break;
            }
            case TokType::word:
            {
                node.kind = Dictionary::Kind::word;
                node.word = t.text;
                expectSemi();
                break;
            }
            case TokType::lbrace:
            {
                node.kind = Dictionary::Kind::block;
                node.block = std::make_shared<Dictionary>(parseBlockBody());
                break;
            }
            case TokType::lparen:
            {
                node.kind = Dictionary::Kind::list;
                parseListBody(node);
                break;
            }
            default:
                cfFatal("dictionary: bad value for key '" + key + "'");
        }
        dict.insert(key, std::move(node));
    }

    void parseListBody(Dictionary::Node& node)
    {
        for (;;)
        {
            Token t = lex_.next();
            if (t.type == TokType::rparen) break;
            if (t.type == TokType::eof) cfFatal("dictionary: unterminated list");
            node.words.push_back(t.text);
        }
        bool allNumeric = !node.words.empty();
        for (const auto& w : node.words)
        {
            char* end = nullptr;
            const double v = std::strtod(w.c_str(), &end);
            if (end != w.c_str() + w.size()) { allNumeric = false; break; }
            node.numbers.push_back(v);
        }
        if (!allNumeric) node.numbers.clear();
        expectSemi();
    }

    // Parses entries until the closing '}' of a block.
    Dictionary parseBlockBody()
    {
        Dictionary block;
        for (;;)
        {
            Token key = lex_.next();
            if (key.type == TokType::eof) cfFatal("dictionary: unterminated block");
            if (key.type == TokType::rbrace) return block;
            if (key.type != TokType::word)
            {
                cfFatal("dictionary: expected block key, got '" + key.text + "'");
            }
            parseEntry(block, key.text);
        }
    }

    void expectSemi()
    {
        Token t = lex_.next();
        if (t.type != TokType::semi)
        {
            cfFatal("dictionary: expected ';' near '" + t.text + "'");
        }
    }

    Lexer lex_;
};

} // namespace

void Dictionary::insert(const std::string& key, Node node)
{
    if (!store_.count(key))
    {
        keys_.push_back(key);
    }
    store_[key] = std::move(node);
}

Dictionary Dictionary::parseString(const std::string& text)
{
    Parser p(text);
    return p.parseTop();
}

Dictionary Dictionary::parseFile(const std::string& path)
{
    std::ifstream in(path);
    if (!in)
    {
        cfFatal("dictionary: cannot open file '" + path + "'");
    }
    std::ostringstream ss;
    ss << in.rdbuf();
    return parseString(ss.str());
}

const Dictionary::Node& Dictionary::find(const std::string& key) const
{
    auto it = store_.find(key);
    if (it == store_.end())
    {
        cfFatal("dictionary: missing required key '" + key + "'");
    }
    return it->second;
}

double Dictionary::getScalar(const std::string& key) const
{
    const Node& n = find(key);
    if (n.kind == Kind::scalar) return n.scalar;
    if (n.kind == Kind::word)
    {
        char* end = nullptr;
        const double v = std::strtod(n.word.c_str(), &end);
        if (end == n.word.c_str() + n.word.size()) return v;
    }
    cfFatal("dictionary: key '" + key + "' is not a scalar");
}

int Dictionary::getInt(const std::string& key) const
{
    return static_cast<int>(getScalar(key));
}

std::string Dictionary::getWord(const std::string& key) const
{
    const Node& n = find(key);
    if (n.kind != Kind::word && n.kind != Kind::scalar)
    {
        cfFatal("dictionary: key '" + key + "' is not a word");
    }
    return n.word;
}

bool Dictionary::getBool(const std::string& key) const
{
    const std::string w = getWord(key);
    if (w == "yes" || w == "true" || w == "on") return true;
    if (w == "no" || w == "false" || w == "off") return false;
    cfFatal("dictionary: key '" + key + "' is not boolean");
}

std::vector<double> Dictionary::getScalarList(const std::string& key) const
{
    const Node& n = find(key);
    if (n.kind != Kind::list || n.numbers.empty())
    {
        cfFatal("dictionary: key '" + key + "' is not a numeric list");
    }
    return n.numbers;
}

std::vector<std::string> Dictionary::getWordList(const std::string& key) const
{
    const Node& n = find(key);
    if (n.kind != Kind::list)
    {
        cfFatal("dictionary: key '" + key + "' is not a list");
    }
    return n.words;
}

const Dictionary& Dictionary::getDict(const std::string& key) const
{
    const Node& n = find(key);
    if (n.kind != Kind::block || !n.block)
    {
        cfFatal("dictionary: key '" + key + "' is not a block");
    }
    return *n.block;
}

double Dictionary::getScalarOr(const std::string& key, double def) const
{
    return has(key) ? getScalar(key) : def;
}

std::string Dictionary::getWordOr(const std::string& key, const std::string& def) const
{
    return has(key) ? getWord(key) : def;
}

} // namespace cf
