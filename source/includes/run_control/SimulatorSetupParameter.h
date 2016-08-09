#ifndef SIMULATORSETUPPARAMETER_H
#define SIMULATORSETUPPARAMETER_H

#include "Parameter.h"

namespace csmp {

/**
    Relates to the geothermal simulator application
    
    @todo documentation is needed (probably Julian Mindel)
    
    @todo SKM move to applications, geothermal energy storage simulator
*/
struct SimulatorSetupParameter : Parameter {
    SimulatorSetupParameter()
    {
        name="Undefined";
        notation="Undefined";
        unit="Undefined";
        key=csmp::Index();
        min=std::numeric_limits<double64>::min();
        max=std::numeric_limits<double64>::max();
        usage="Undefined";
        explanation="Undefined";
        type=SCALAR;
        vsize=1;
        placement=UNDEFINED;
        reference="Undefined";
    }
    SimulatorSetupParameter( const SimulatorSetupParameter& p){
        this->name        = p.name;
        this->notation    = p.notation;
        this->unit        = p.unit;
        this->key         = p.key;
        this->min         = p.min;
        this->max         = p.max;
        this->usage       = p.usage;
        this->explanation = p.explanation;
        this->type        = p.type;
        this->vsize       = p.vsize;
        this->placement   = p.placement;
        this->reference   = p.reference;

    }
    SimulatorSetupParameter& operator=( const SimulatorSetupParameter& p)
    {
        if ( &p != this ) {
            this->name        = p.name;
            this->notation    = p.notation;
            this->unit        = p.unit;
            this->key         = p.key;
            this->min         = p.min;
            this->max         = p.max;
            this->usage       = p.usage;
            this->explanation = p.explanation;
            this->type        = p.type;
            this->vsize       = p.vsize;
            this->placement   = p.placement;
            this->reference   = p.reference;
        }
        return *this;
    }
    bool operator==( const SimulatorSetupParameter& p ) const
    {
        return (p.name==name &&
                p.notation==notation &&
                p.unit==unit &&
                p.type==type &&
                p.vsize==p.vsize &&
                p.min==min &&
                p.max==max &&
                p.placement==placement &&
                p.usage==usage);
    }

    bool operator!=( const SimulatorSetupParameter& p ) const
    {
        return (p.name!=name ||
                p.notation!=notation ||
                p.unit!=unit ||
                p.type!=type ||
                p.vsize!=p.vsize ||
                p.min!=min ||
                p.max!=max ||
                p.placement!=placement ||
                p.usage!=usage);
    }

    ~SimulatorSetupParameter(){}

    VARIABLE_TYPE    type;
    size_t           vsize;
    PLACEMENT        placement;
};

bool compare_setup_parameter_nocase (const csmp::SimulatorSetupParameter& first, const csmp::SimulatorSetupParameter& second);

} // end csmp

#endif // SIMULATORSETUPPARAMETER_H

