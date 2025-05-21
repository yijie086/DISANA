#include "Event.h"
#include <iostream>

Event::Event(const std::string& fileName)
    : clas12::clas12reader(fileName), validEvent(false) { 
    std::cout << "Opening file: " << getFilename() << std::endl;
    std::cout << "Is file open? " << (this->isOpen() ? "Yes" : "No") << std::endl;

    if (!this->isOpen()) {
        std::cerr << "Error: Failed to open file: " << fileName << std::endl;
    }
}

Event::~Event() {
    // No explicit cleanup needed; clas12reader handles file closing internally
}

bool Event::getNextEvent() {
    if (!this->isOpen()) {
        std::cerr << "Error: File is not open!" << std::endl;
        return false;
    }

    validEvent = next();  // Calls clas12reader method to get the next event
    if (validEvent) {
        evt = this->getEvent();  // Properly assign the event
    } else {
        std::cerr << "Warning: No more events available" << std::endl;
    }

    return validEvent;
}

bool Event::isValid() const {
    return validEvent;
}

hipo::event& Event::getEvent() {
    return evt;
}
