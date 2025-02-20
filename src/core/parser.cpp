#include "parser.h"

#include <iostream>

#include <boost/algorithm/string.hpp>

using std::cout;
using std::endl;
using std::string;
using std::vector;

Parser::Parser(const string filename, const FileFormat format)
{
    // TODO preprocess file to include ITPs
    format_   = format;
    filename_ = filename;
    file_.open(filename);
    if (!file_.is_open())
        throw std::runtime_error("File " + filename + " could not be opened");
}

Parser::~Parser()
{
    if (file_.is_open())
        file_.close();
}

bool Parser::getLine(vector<string> &tokens)
{
    while (true)
    {
        eof_ = !getline(file_, line_);

        // Stop if we hit eof
        if (eof_)
            return false;
        boost::trim(line_);

        // Skip comments
        if (line_[0] == ';' || line_[0] == '#')
            continue;

        // Line is empty, ignore it
        if (line_ == "")
            continue;

        switch (format_)
        {
            case FileFormat::GROMACS:
                // Line is a section header
                if (line_[0] == '[')
                {
                    section_ = line_.substr(line_.find_first_of('[') + 1,
                                            line_.find_last_of(']') - 1);
                    boost::trim(section_);
                    continue;
                }
                break;

            case FileFormat::LAMMPS:
                throw std::logic_error("Not implemented");
        };

        // Line isn't empty, accept it
        break;
    }

    // Separate and trim whitespace from tokens
    vector<string> tmp;
    boost::split(tmp, line_, boost::is_any_of(";"),
                 boost::algorithm::token_compress_on);
    boost::trim(tmp[0]);
    boost::split(tokens, tmp[0], boost::is_any_of("\t "),
                 boost::algorithm::token_compress_on);
    for (string &tok : tokens)
        boost::trim(tok);

    // Return true if there is still file to read
    return true;
}

bool Parser::findSection(const string find)
{
    rewind();
    vector<string> token_buffer;
    while (section_ != find)
    {
        if (!getLine(token_buffer))
            return false;
    }
    rewind();
    return true;
}

bool Parser::getLineFromSection(const string find, vector<string> &tokens,
                                const int len)
{
    // Are we looking for a new section? - it might be above the last one
    if (find != findPrevious_)
        rewind();
    findPrevious_ = find;
    while (getLine(tokens))
    {
        if (section_ == find && tokens.size() >= len)
            return true;
    }
    rewind();
    return false;
}

void Parser::rewind()
{
    // Clear eof and rewind
    eof_ = false;
    file_.clear();
    file_.seekg(0, std::ios::beg);
}

bool Parser::getKeyFromSection(const string &section, const string &key,
                               string &value)
{
    vector<string> tmp;
    while (getLineFromSection(section, tmp, 2))
    {
        if (tmp[0] == key)
        {
            value = tmp[1];
            return true;
        }
    }
    return false;
}

int Parser::getIntKeyFromSection(const string &section, const string &key,
                                 const int default_value)
{
    rewind();
    string tmp;
    if (getKeyFromSection(section, key, tmp))
        return stoi(tmp);
    return default_value;
}

double Parser::getDoubleKeyFromSection(const string &section, const string &key,
                                       const double default_value)
{
    rewind();
    string tmp;
    if (getKeyFromSection(section, key, tmp))
        return stof(tmp);
    return default_value;
}

string Parser::getStringKeyFromSection(const string &section, const string &key,
                                       const string &default_value)
{
    rewind();
    string tmp;
    if (getKeyFromSection(section, key, tmp))
    {
        boost::to_upper(tmp);
        return tmp;
    }
    return default_value;
}
