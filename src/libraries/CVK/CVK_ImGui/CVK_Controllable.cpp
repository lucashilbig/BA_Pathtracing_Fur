#include "CVK_Controllable.h"

CVK::Controllable::Controllable(std::string title)
{
    m_title = title;
}

CVK::Controllable::~Controllable()
{

}

void CVK::Controllable::updateGui()
{
    if(ImGui::CollapsingHeader(m_title.c_str()))
    {
        fillGui();
    }
}
