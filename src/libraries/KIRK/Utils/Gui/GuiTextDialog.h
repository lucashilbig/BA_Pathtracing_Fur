#ifndef __KIRK_GUI_TEXT_DIALOG
#define __KIRK_GUI_TEXT_DIALOG

#include <string>
#include "Gui.h"

namespace KIRK
{
	class GuiTextDialog : public GuiNamedElement
	{
	public:
		GuiTextDialog(const std::string &title, const std::string &message);

		void onGui() override;

		GuiTextDialog &addPositive(std::string positive = "OK", std::function<void()> onclick = nullptr);
		GuiTextDialog &addNegative(std::string negative = "Cancel", std::function<void()> onclick = nullptr);

	private:
		std::function<void()> m_onclick_positive;
		std::function<void()> m_onclick_negative;

		bool m_has_positive = false;
		bool m_has_negative = false;

		bool m_opened = false;

		std::string m_content;
		std::string m_positive;
		std::string m_negative;
	};

}

#endif // !__KIRK_GUI_TEXT_DIALOG
