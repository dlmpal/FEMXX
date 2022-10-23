#ifndef ELEMENTMANAGER_HPP 
#define ELEMENTMANAGER_HPP 


#include "Element.hpp"
#include "Elements.hpp"
#include <vector>


class ElementManager{

public:
    ElementManager(){};
    ElementManager(int NGP , MeshData& meshData);
    ElementManager(const ElementManager& rhs);
    ElementManager& operator=(const ElementManager& rhs);
    ~ElementManager();

private: 
    void BuildElements(MeshData& meshData);

private: 
    int NumEle; 
    int NGP; // Number of Gauss Points
    std::vector<Element> FE; 
    ElementPrototype FE_PROTO; 
    friend class Model;
};



#endif // ELEMENTMANAGER_HPP