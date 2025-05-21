#include "EventCuts.h"
#include <iostream>

EventCuts::EventCuts() : validEvent(false) {}
EventCuts::EventCuts(const EventCuts &cuts) {}

EventCuts &EventCuts::operator=(const EventCuts &cuts) {
  if (this == &cuts) {
    return *this;
  }
  return *this;
}
EventCuts::~EventCuts() {}
bool EventCuts::isSelected(Event *evt) {
  bool pass = true;
  /// here I will define the creteria for the selection of some particles
  /// related to processes
  // evt->SetSelectionStatus(pass);
  return pass;
}