#ifndef __KIRK_GUI_SCENE_H
#define __KIRK_GUI_SCENE_H



#include "Gui.h"
#include <KIRK/Common/SceneGraph.h>
#include "../../../CVK/CVK_Utils/CVK_ConverterUtils.h"

namespace KIRK {
	enum SceneUpdateFlags
	{
		UPDATE_MATERIALS = 1 << 0,
		UPDATE_GRAPH = 1 << 1
	};

class GuiScene
{
public:
	GuiScene(std::shared_ptr<KIRK::Gui> gui);
	~GuiScene();

	void buildSceneGui(std::shared_ptr<SceneGraph> graph, CVK::CVKCameraSynchronizer &sync);
	void setUpdateCallback(std::function<void(unsigned)> callback, unsigned not_editable = 0);

private:
	std::shared_ptr<KIRK::Gui> m_gui;
	std::function<void(unsigned)> m_update_callback;
	std::vector<std::shared_ptr<Camera>> m_cameras;
	int m_current_camera = 0;
	unsigned m_flags;
	unsigned m_editable_flags;
};
}

#endif // !__KIRK_GUI_SCENE_H
