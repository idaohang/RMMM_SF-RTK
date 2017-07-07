//
// Created by sujinglan on 7/7/17.
//

/**
 * @file arc_MovementModel.cc
 * @brief Source file for the ARC_MovementModel class.
 * @author sujinglan
 */

#include "arc_MovementModel.h"
#include "arc_PF.h"
#include <libPF/CRandomNumberGenerator.h>

ARC_MovementModel::ARC_MovementModel() : libPF::MovementModel<ARC_States>() ,
                                         m_Std_BaseX(ARC_PF_BASEPOS_STD),
                                         m_Std_BaseY(ARC_PF_BASEPOS_STD),
                                         m_Std_BaseZ(ARC_PF_BASEPOS_STD),
                                         m_Std_RoverX(ARC_PF_ROVERPOS_STD),
                                         m_Std_RoverY(ARC_PF_ROVERPOS_STD),
                                         m_Std_RoverZ(ARC_PF_ROVERPOS_STD),
                                         m_Std_Rover_Clk(ARC_PF_ROVERCLK_STD),
                                         m_Std_Base_Clk(ARC_PF_BASECLK_STD),
                                         m_Std_Base_Trp(ARC_PF_BASETROP_STD),
                                         m_Std_Rover_Trp(ARC_PF_ROVERTROP_STD)
{
    m_RNG = new libPF::CRandomNumberGenerator();
    for (size_t i=0;i<MAXSAT;i++) m_Std_N[i]=ARC_PF_AMB_STD;
}

ARC_MovementModel::~ARC_MovementModel() {
    if (m_RNG) delete m_RNG;
}