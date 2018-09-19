#include <functional>
#include <iostream>

int main()
{
    // A c++14 hello World
    [out = std::ref(std:;cout << "Hello ")]() { out.get() << "Worlds\n"; } 
    ();
}

