//******************************************************************************
// Config.h
//	: header file for system configuration
//
//------------------------------------------------------------------------------
//
//	Ver 1.00		Date: Tue Jun 19 02:36:37 2019
//
//******************************************************************************

#ifndef _Base_Config_H
#define _Base_Config_H

//------------------------------------------------------------------------------
//	Including Header Files
//------------------------------------------------------------------------------

#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
#ifndef _WIN32
    #include <unistd.h>
#endif


using namespace::std;

//------------------------------------------------------------------------------
//	Defining Macros
//------------------------------------------------------------------------------

namespace ZuKai {
namespace Base {

    //------------------------------------------------------------------------------
    //	Defining Classes
    //------------------------------------------------------------------------------
	class Config {

    private:

        map< string, string > m_map;

        //------------------------------------------------------------------------------
        //	Special functions
        //------------------------------------------------------------------------------
        void split( const string& in, vector< string >& out );
        void loadConfigFile( const string& filename );

    protected:

    public:

        //------------------------------------------------------------------------------
        //	Constructors & Destructors
        //------------------------------------------------------------------------------
        // default constructor
        Config( void );
        // parameterized constructor
        Config( const string& filename = "" ) {
            loadConfigFile( filename );
        }
        // copy constructor
        Config( const Config & c ) {}
        // destructor
        ~Config( void ){}

        //------------------------------------------------------------------------------
        //	Special functions
        //------------------------------------------------------------------------------
        double      getlf   ( const string& key ) const;
        float       getf    ( const string& key ) const;
        int         geti    ( const string& key ) const;
        string      gets    ( const string& key ) const;
        bool        has     ( const string& key ) const;

        //------------------------------------------------------------------------------
        //	I/O functions
        //------------------------------------------------------------------------------
        // output
        friend ostream &	operator << ( ostream & s, const Config & c );
        // input
        friend istream &	operator >> ( istream & s, Config & c );
        // class name
        virtual const char * className( void ) const { return "Config"; }
    };

} // namespace Base
} // namespace ZuKai


#endif // _Base_Config_H
