#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "cf/common/Error.h"

namespace cf
{

/// Lightweight nested key/value text database, inspired by OpenFOAM's
/// dictionary class but with an independent, simpler grammar:
///
///   key 3.14;
///   name myWall;
///   text "quoted string";
///   list (1 2 3);
///   words (inlet outlet);
///   block
///   {
///       type fixedValue;
///       value 0;
///   }
///
/// '#' and '//' start line comments.
class Dictionary
{
public:
    enum class Kind { scalar, word, list, block };

    struct Node
    {
        Kind kind = Kind::word;
        double scalar = 0;
        std::string word;
        std::vector<double> numbers;
        std::vector<std::string> words;
        std::shared_ptr<Dictionary> block;
    };

    Dictionary() = default;
    Dictionary(Dictionary&&) = default;
    Dictionary& operator=(Dictionary&&) = default;
    Dictionary(const Dictionary&) = default;
    Dictionary& operator=(const Dictionary&) = default;

    static Dictionary parseFile(const std::string& path);
    static Dictionary parseString(const std::string& text);

    bool has(const std::string& key) const { return store_.count(key) > 0; }
    std::size_t size() const { return keys_.size(); }
    const std::vector<std::string>& keys() const { return keys_; }

    double getScalar(const std::string& key) const;
    int getInt(const std::string& key) const;
    std::string getWord(const std::string& key) const;
    bool getBool(const std::string& key) const;
    std::vector<double> getScalarList(const std::string& key) const;
    std::vector<std::string> getWordList(const std::string& key) const;
    const Dictionary& getDict(const std::string& key) const;
    const Node& find(const std::string& key) const;

    double getScalarOr(const std::string& key, double def) const;
    std::string getWordOr(const std::string& key, const std::string& def) const;

    void insert(const std::string& key, Node node);

private:
    std::map<std::string, Node> store_;
    std::vector<std::string> keys_; // distinct keys, first-occurrence order
};

} // namespace cf
