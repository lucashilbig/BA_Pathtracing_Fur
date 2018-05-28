#ifndef CVK_CONTROLLABLE_H
#define CVK_CONTROLLABLE_H

#include <string>
#include <externals/ImGui/imgui.h>

namespace CVK {
class Controllable
{
public:
    Controllable(std::string title);
    virtual ~Controllable();
    void updateGui();
protected:
    virtual void fillGui() {}
	std::string m_title;
};
}
#endif // CVK_CONTROLLABLE_H
