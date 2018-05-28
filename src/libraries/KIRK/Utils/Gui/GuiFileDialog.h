#ifndef __KIRK_GUI_FILEDIALOG_H
#define __KIRK_GUI_FILEDIALOG_H

#include <string>
#include <experimental/filesystem>
#include <memory>
#include <functional>

#include "Gui.h"

namespace fs = std::experimental::filesystem;

namespace KIRK
{

	class GuiFileDialog : public GuiNamedElement, public std::enable_shared_from_this<GuiFileDialog>
	{
	public:
		using dialog_submit_callback = std::function<void(const GuiFileDialog &dialog, const fs::path &path)>;
		GuiFileDialog(const std::string &title, const std::string &positive_action, const std::string &negative_action, const std::string &start_path, bool check_if_exists, dialog_submit_callback on_submit);
		~GuiFileDialog() override;

		void onGui() override;

	private:
		bool m_opened = false;
		fs::path m_current_path;
		std::string m_filename;
		bool m_check_if_exists;
		bool m_file_exists = false;
		std::string m_title;
		std::string m_positive;
		std::string m_negative;
		dialog_submit_callback m_submit_callback;
	};
}

#endif // !__KIRK_GUI_FILEDIALOG_H
