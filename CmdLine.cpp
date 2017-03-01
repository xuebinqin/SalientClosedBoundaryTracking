#include <cstdlib>
#include <cstring>
using namespace std;

#include "CmdLine.h"

CmdLine::CmdLine(int argc, char *argv[])
    : _moreOptions(true)
{
    /*
    // Redisplay the command line, as received from the shell
    for(int i = 0; i < argc; i++)
    {
	cout << argv[i] << " ";
    }
    cout << endl;
    */

    // Find the "--" on the command line (if any)
    indexOfDashDash = 1;
    while(indexOfDashDash < argc && string("--") != argv[indexOfDashDash])
    {
	indexOfDashDash++;
    }

    // Store all of the command line arguments into a vector of strings
    {
	int i;
	for(i = 1; i < indexOfDashDash; i++)
	{
	    int len = strlen(argv[i]);
	    if (len > 1 && argv[i][0] == '-' && argv[i][1] != '-')
	    {
		for (int j = 1; j < len; j++)
		{
		    args.push_back(string("-") + argv[i][j]);
		}
	    }
	    else
	    {
		args.push_back(argv[i]);
	    }
	}
	indexOfDashDash = args.size();
	for(; i < argc; i++)
	{
	    args.push_back(argv[i]);
	}
    }

    /*
    // Display the arguments, as parsed so far
    cout << indexOfDashDash << endl;
    for(int i = 0; i < args.size(); i++)
    {
	cout << args[i] << endl;
    }
    */
}

CmdArg CmdLine::parse(char shortoption, const string &longoption)
{
    bool found = false;

    string foundStr = extract(shortoption, longoption, false);
    if(foundStr != "")
    {
	found = true;
    }
    return CmdArg(found);
}

CmdArgInt CmdLine::parse(char shortoption, const string &longoption, int defaultValue)
{
    bool found = false;
    int value = defaultValue;

    string valueStr = extract(shortoption, longoption, true);
    if(valueStr != "")
    {
	found = true;
	//--- There needs to be some type checking here
	value = atoi(valueStr.c_str());
    }
    return CmdArgInt(found, value);
}

CmdArgDouble CmdLine::parse(char shortoption, const string &longoption, double defaultValue)
{
    bool found = false;
    double value = defaultValue;

    string valueStr = extract(shortoption, longoption, true);
    if(valueStr != "")
    {
	found = true;
	//--- There needs to be some type checking here
	value = atof(valueStr.c_str());
    }
    return CmdArgDouble(found, value);
}

CmdArgString CmdLine::parse(char shortoption, const string &longoption, const string &defaultValue)
{
    bool found = false;
    string value = defaultValue;

    string valueStr = extract(shortoption, longoption, true);
    if(valueStr != "")
    {
	found = true;
	value = valueStr;
    }
    return CmdArgString(found, value);
}

string CmdLine::extract(char shortoption, const string &longoption, bool takesValue)
{
    string value("");                                // null string means no match was found for option
    for(int i = 0; i < indexOfDashDash && value == ""; i++)
    {
	if(args[i].size() > 1 && args[i][0] == '-')
	{
	    if(args[i][1] == '-')                        // it is a long option
	    {
		// Find the "=" (if any)
		int indexOfEquals = 2;
		while(indexOfEquals < (int)args[i].size() && args[i][indexOfEquals] != '=')
		{
		    indexOfEquals++;
		}
		if(indexOfEquals > 2)
		{
		    if(longoption == args[i].substr(2, indexOfEquals - 2))  // In other words, if they are a match
		    {
			if(takesValue)
			{
			    if(indexOfEquals + 1 < (int)args[i].size())
			    {
				value = args[i].substr(indexOfEquals + 1, args[i].size() - indexOfEquals - 1);
				args[i] = "";
			    }
			    else
			    {
				//--- There needs to be a less sloppy errorMessage
				string errorMessage;
				errorMessage = errorMessage + "Option " + shortoption + " " + longoption + " requires a value: " + args[i];
				throw errorMessage;
			    }
			}
			else
			{
			    if(indexOfEquals == (int)args[i].size())
			    {
				value = "true";
				args[i] = "";
			    }
			    else
			    {
				//--- There needs to be a less sloppy errorMessage
				string errorMessage;
				errorMessage = errorMessage + "Option " + shortoption + " " + longoption + " does not take a value: " + args[i];
				throw errorMessage;
			    }
			}
		    }
		}
	    }
	    else                                         // it is a shortoption
	    {
		if(args[i][1] == shortoption)
		{
		    if(takesValue)
		    {
			if(i + 1 < indexOfDashDash && args[i + 1][0] != '-')
			{
			    value = args[i + 1];
			    args[i] = "";
			    args[i + 1] = "";
			}
			else
			{
			    //--- There needs to be a less sloppy errorMessage
			    string errorMessage;
			    errorMessage = errorMessage + "Option " + shortoption + " " + longoption + " requires a value: " + args[i];
			    throw errorMessage;
			}
		    }
		    else
		    {
			value = "true";
			args[i] = "";
		    }
		}
	    }
	}
    }

    /*
    // Display a confirmation message
    cout << shortoption << "," << longoption << " " << (value != "" ? "found..." : "not found...");
    cout << "value: " << value << endl;
    */

    return value;
}

CmdArgString CmdLine::getNextOperand()
{
    if(moreOptions())
    {
	string errorMessage;
	errorMessage = errorMessage + "Leftover options on command line. Cannot process request for next operand.";
	throw errorMessage;
    }

    bool found = false;
    string value("");

    int i = 0;
    while(i < (int)args.size() && (args[i] == "" || args[i] == "--"))
    {
	i++;
    }
    if(i < (int)args.size())
    {
	found = true;
	value = args[i];
	args[i] = "";
    }

    /*
    // Display a confirmation message
    cout << "operand " << (value != "" ? "found..." : "not found...");
    cout << "value: " << value << endl;
    */

    return CmdArgString(found, value);
}

bool CmdLine::moreOptions()
{
    if(_moreOptions)          // We only need to check if we have not already determined that there are no more options
    {
	// Let's double check
	_moreOptions = false;
	for(int i = 0; i < indexOfDashDash && !_moreOptions; i++)
	{
	    if(args[i].size() > 1 && args[i][0] == '-')
	    {
		_moreOptions = true;
	    }
	}
    }

    return _moreOptions;
}

