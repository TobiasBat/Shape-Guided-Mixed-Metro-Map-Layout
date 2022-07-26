//******************************************************************************
// Config.cpp
//	: program file for system configuration
//
//------------------------------------------------------------------------------
//
//	Ver 1.00		Date: Sun Aug 16 15:02:45 2020
//
//******************************************************************************

//------------------------------------------------------------------------------
//	Including Header Files
//------------------------------------------------------------------------------

#include "Common.h"

using namespace std;

namespace ZuKai {
namespace Base {
	
	double Common::_mainwidget_width = 800;
	double Common::_mainwidget_height = 600;
	double Common::_dockwidget_width = 60;
	double Common::_menubar_height = 60;
	
	//------------------------------------------------------------------------------
	//	Private Functions
	//------------------------------------------------------------------------------
	
	//------------------------------------------------------------------------------
	//	Protected Functions
	//------------------------------------------------------------------------------
	
	//------------------------------------------------------------------------------
	//	Public functions
	//------------------------------------------------------------------------------
	//  Common::_stringToDouble -- convert string to double
	//
	//  Inputs
	//      string
	//
	//  Outputs
	//  double
	//
	double Common::stringToDouble( const std::string str )
	{
		stringstream ss( str );
		double val = 0;
		ss >> val;
		
		return val;
	}
	
	
	const double &Common::getMainwidgetWidth() {
		return _mainwidget_width;
	}
	
	void Common::setMainwidgetWidth( double mainwidgetWidth ) {
		_mainwidget_width = mainwidgetWidth;
	}
	
	const double &Common::getMainwidgetHeight() {
		return _mainwidget_height;
	}
	
	void Common::setMainwidgetHeight( double mainwidgetHeight ) {
		_mainwidget_height = mainwidgetHeight;
	}
	
	const double Common::getDockWidgetWidth() {
		return _dockwidget_width;
	}
	
	void Common::setDockWidgetWidth( double dockwidgetWidth ) {
		_dockwidget_width = dockwidgetWidth;
	}
	
	const double Common::getMenubarHeight() {
		return _menubar_height;
	}
	
	void Common::setMenubarHeight( double menubarHeight ) {
		_menubar_height = menubarHeight;
	}

} // namespace Base
} // namespace ZuKai
