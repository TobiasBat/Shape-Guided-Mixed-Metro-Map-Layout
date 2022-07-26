//******************************************************************************
// Config.cpp
//	: program file for system configuration
//
//------------------------------------------------------------------------------
//
//	Ver 1.00		Date: Sun Sep 16 15:02:45 2019
//
//******************************************************************************

//------------------------------------------------------------------------------
//	Including Header Files
//------------------------------------------------------------------------------

#include "Config.h"

using namespace std;

namespace ZuKai {
namespace Base {

    //------------------------------------------------------------------------------
    //	Private Functions
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    //	Protected Functions
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    //	Public functions
    //------------------------------------------------------------------------------
    //
    //  Config::getlf --	string to double
    //
    //  Inputs
    //	key	: text string
    //
    //  Outputs
    //	the corresponding value
    //
    double Config::getlf( const string& key ) const {

        map< string, string >::const_iterator it = m_map.find( key );

        if( it == m_map.end() )
        {
            cerr << "[ERR] " << key << " is not a defined option." << endl;
            return 0.0;
        }

        return atof( it->second.c_str() );
    }

    //
    //  Config::getf --	string to float
    //
    //  Inputs
    //	key	: text string
    //
    //  Outputs
    //	the corresponding value
    //
    float Config::getf( const string& key ) const
    {
        map< string, string >::const_iterator it = m_map.find( key );

        if( it == m_map.end() )
        {
            cerr << "[ERR] " << key << " is not a defined option." << endl;
            return 0.0f;
        }

        return atof( it->second.c_str() );
    }

    //
    //  Config::geti --	string to integer
    //
    //  Inputs
    //	key	: text string
    //
    //  Outputs
    //	the corresponding value
    //
    int Config::geti( const string& key ) const
    {
        map< string, string >::const_iterator it = m_map.find(key);

        if( it == m_map.end() )
        {
            cerr << "[ERR] " << key << " is not a defined option." << endl;
            return 0;
        }

        return atoi( it->second.c_str() );
    }

    //
    //  Config::gets --	string to string
    //
    //  Inputs
    //	key	: text string
    //
    //  Outputs
    //	the corresponding value
    //
    string Config::gets( const string& key ) const
    {
        map< string, string>::const_iterator it = m_map.find( key );

        if( it == m_map.end() )
        {
            cerr << "[ERR] " << key << " is not a defined option." << endl;
            return string( "" );
        }

        return it->second;
    }

    //
    //  Config::has --	check if in the list
    //
    //  Inputs
    //	key	: text string
    //
    //  Outputs
    //	yes/no
    //
    bool Config::has( const string& key ) const
    {
        if( m_map.find( key ) == m_map.end() )
        {
            return false;
        }

        return true;
    }

    //
    //  Config::split --	separate by " "
    //
    //  Inputs
    //	in	: text string
    //	out	: text string vector
    //
    //  Outputs
    //	none
    //
    void Config::split( const string& in, vector< string >& out )
    {
        out.clear();
        out.resize(2);

        string::size_type sp_st  = in.find( " ", 0 );
        string::size_type sp_ed  = sp_st;

        if( sp_st == string::npos ) return;

        while( true ){
            if( in.at( ++sp_ed ) != ' ' )
            {
                break;
            }
        }

        out[0] = in.substr( 1, sp_st - 1 );
        out[1] = in.substr( sp_ed, in.length() );

    }

    //
    //  Config::loadConfigFile --	load file
    //
    //  Inputs
    //	filename	: text string
    //
    //  Outputs
    //	none
    //
    void Config::loadConfigFile( const string& filename )
    {
        ifstream ifs( filename.c_str(), ios::in );
        if( ifs.fail() )
        {
            cerr << "[ERR] file open error " << filename << endl;
            return;
        }

        string buf;

        while( getline( ifs, buf ) )
        {
            if( buf.length() == 0   ) continue;
            if( buf.at(0)    != ':' ) continue;

            vector< string > vbuf;
            split( buf, vbuf );
    #ifdef  DEBUG
            cerr << "buf = " << buf << endl;
                        for( int i = 0; i < vbuf.size(); i++ ){
                            cerr << "vbuf[" << i << "] = " << vbuf[i] << endl;
                        }
    #endif  // DEBUG
            m_map.insert( make_pair(vbuf[0], vbuf[1]) );
        }
    }

} // namespace Base
} // namespace ZuKai
