//
// Created by sujinglan on 7/7/17.
//

#ifndef ARC_ARC_MOVEMENTMODEL_H
#define ARC_ARC_MOVEMENTMODEL_H
#include <libPF/MovementModel.h>
#include "arc_States.h"

/**
 * @class ARC-SRTK MovementModel
 *
 * @brief single frequency rtk position movement model propagates class
 *
 * This movement model propagates a car state according to the translation and
 * rotation speed.
 *
 * @author sujinglan
 */
class ARC_MovementModel : public libPF::MovementModel<ARC_States>
{
    /**
     * Constructor
     */
    ARC_MovementModel();

    /**
     * Destructor
     */
    ~ARC_MovementModel();

    /**
     * The drift method propagates the car using its speed.
     * @param state Pointer to the state that has to be manipulated.
     */
    void drift(ARC_States& state, double dt) const;

    /**
     * The diffusion consists of a very small gaussian jitter on the
     * state's variable.
     * @param state Pointer to the state that has to be manipulated.
     */
    void diffuse(ARC_States& state, double dt) const;

    /**
     * @param d new standard deviation for the diffusion of x
     */
    void setXStdDev(double d);

private:
    
};
#endif //ARC_ARC_MOVEMENTMODEL_H
