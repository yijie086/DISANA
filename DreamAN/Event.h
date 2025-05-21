#ifndef EVENT_H
#define EVENT_H

#include "hipo4/event.h"
#include "clas12reader.h"
#include "particle.h"

class Event : public clas12::clas12reader {
public:
    // Constructor initializes the base class (clas12reader)
    explicit Event(const std::string& fileName);

    // Destructor
    ~Event();

    // Fetches the next event, returning true if successful
    bool getNextEvent();

    // Checks if event retrieval was successful
    bool isValid() const;

    // Gets the current event
    hipo::event& getEvent();

private:
    bool validEvent;  // Flag to check if event retrieval was successful
    hipo::event evt;  // Stores the retrieved event
};

#endif // EVENT_H
