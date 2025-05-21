#ifndef EVENTCUTS_H
#define EVENTCUTS_H

#include "Event.h"
#include "hipo4/event.h"

class EventCuts {
public:
  EventCuts();
  EventCuts(const EventCuts &cuts);
  EventCuts &operator=(const EventCuts &cuts); // Constructor that takes a
                                               // reference to DataReader
  virtual ~EventCuts();
  bool isSelected(Event *evt);

private:
  bool validEvent; // Flag to check if event retrieval was successful
};

#endif // EVENT_H
