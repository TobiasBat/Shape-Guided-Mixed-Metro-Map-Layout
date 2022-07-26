#pragma once
#include <QtWidgets/QWidget>

namespace Ui
{
	class BaseView
	{
	public:
		BaseView();
		virtual ~BaseView();
		virtual QWidget* getSettings() = 0;
		
		void loadFile();

	protected:
		QWidget* m_settings;

		virtual void processFile(std::string path) = 0;
	};
}