#ifndef _CmdLine_h_
#define _CmdLine_h_

#include <vector>
#include <string>
#include <iostream>
#include <string>
using namespace std;

class CmdArg
{
public:
    CmdArg(bool = false);
    bool isFound() const;

private:
    bool _isFound;
};

inline CmdArg::CmdArg(bool found) : _isFound(found) { }
inline bool CmdArg::isFound() const { return _isFound; }

class CmdArgInt : public CmdArg
{
public:
    CmdArgInt(bool = false, int = -999);
    int value() const;
    void setValue(int);

private:
    int _value;
};

inline CmdArgInt::CmdArgInt(bool found, int defaultValue) : CmdArg(found), _value(defaultValue) { }
inline int CmdArgInt::value() const { return _value; }
inline void CmdArgInt::setValue(int i) { _value = i; }

class CmdArgDouble : public CmdArg
{
public:
    CmdArgDouble(bool = false, double = 1e-999);
    double value() const;
    void setValue(double);

private:
    double _value;
};

inline CmdArgDouble::CmdArgDouble(bool found, double defaultValue) : CmdArg(found), _value(defaultValue) { }
inline double CmdArgDouble::value() const { return _value; }
inline void CmdArgDouble::setValue(double d) { _value = d; }

class CmdArgString : public CmdArg
{
public:
    CmdArgString(bool = false, const string & = "");
    const string &value() const;

private:
    string _value;
};

inline CmdArgString::CmdArgString(bool found, const string &defaultValue) : CmdArg(found), _value(defaultValue) { }
inline const string &CmdArgString::value() const { return _value; }

class CmdLine
{
public:
    CmdLine(int, char *[]);

    CmdArg parse(char, const string &);
    CmdArgInt parse(char, const string &, int);
    CmdArgDouble parse(char, const string &, double);
    CmdArgString parse(char, const string &, const string &);
    CmdArgString getNextOperand();
    bool moreOptions();

private:
    bool _moreOptions;
    int indexOfDashDash;
    vector<string> args;

    string extract(char, const string &, bool);
};

#endif // _CmdLine_h_
