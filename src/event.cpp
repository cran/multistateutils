#include "event.h"
#include "simulation.h"
#include "state.h"

#include <tuple>
#include <iostream>

// Constructor
Event::Event(int individual_id, int state, double time, double sim):
    individual_id(individual_id), state_entering(state), time(time), sim_entry(sim) {}

void Event::processEvent(Simulation* sim) {
    sim->add_history(std::tuple<int, int, double> (individual_id, state_entering, time));

    State* entering_state = sim->get_state(state_entering);
    if (entering_state->is_transient()) {

        // Get next transition, as a pair of <state>,<entry time>
        std::pair<int, double> next_transition = entering_state->get_next_transition(sim->get_patient_attrs(individual_id),
                                                                                     (this->time - this->sim_entry));  // Time since entry
        
        sim->add_event(Event(individual_id,
                             next_transition.first,  // state
                             next_transition.second + this->time,  // Transition time is relative, add it to clock
                             this->sim_entry));
    }

}
