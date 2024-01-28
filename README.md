# Nonlinear-ID-and-NMPC-of-the-Lorenz-System
In this file, nonlinear system identification using SINDy(Sparse Identification of Nonlinear Dynamics) is performed to identify the dynamics of the Lorenz System.
Once the system is identified, an NMPC(Nonlinear Model Predictive Control) is applied to predict and control the system to achieve reference tracking.

In this file, you may change input & output references to fit your needs and change the value of threshold when performing STLS (Sequential Total Least Square) in SINDy. Furthermore, you may extend or maunally change the candidate functions in Candidate_Library.m file for your identification task.
