#include "BaseView.h"

namespace Ui
{
	BaseView::BaseView() 
	{
		m_settings = nullptr;
	}
	
	BaseView::~BaseView() 
	{
	}

	void BaseView::loadFile()
	{
		QUrl url = QFileDialog::getOpenFileUrl();

		processFile(url.toLocalFile().toStdString());
	}
}