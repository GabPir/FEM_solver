
# FEM_solver_(from_scratch)

## Project for "Numerical methods for PDEs 2024"

This work is an application of the Galerkin discretization method for both stationary and evolutionary problems. The general model considered is the Reaction-Diffusion-Convection model, initially in its stationary form.
<br>
![immagine](https://github.com/user-attachments/assets/7148bfa8-0a74-410b-8727-b28cbb29cc41)
<br>
The Galerkin method provides an approximate solution to the problem in its variational form, defined in the finite space V_h. The choice of V_h depends on the type of finite element selected and the regularity of the solution. Under appropriate assumptions, the Lax-Milgram theorem guarantees the existence and uniqueness of the approximate solution, as well as its continuous dependence on the data. 
<br>
The extension to the time-dependent case concerns the problem:
<br>
![immagine](https://github.com/user-attachments/assets/0fab8b04-2786-4819-84fa-470ca72188ff)
<br>
and will be addressed using the Method of Lines (MOL), first applying a semi-discretization in space to obtain a system of ODEs in time, which can then be solved using time-stepping discretization methods.


## Stationary Problem Case

As a numerical example, consider the "bubble function" in the domain Ω=[0,1]×[0,1]:

<br>

![immagine](https://github.com/user-attachments/assets/757cba81-8487-4511-8fc0-045d4f846c33)
<br>
<br>

![immagine](https://github.com/user-attachments/assets/31379c71-1359-4287-9a5b-68a9058e4911)
<br>
As for the study of convergence with respect to the maximum mesh diameter, the results of the a priori analysis are confirmed:<br>
![immagine](https://github.com/user-attachments/assets/97bb3693-614e-4152-abb4-143542d58540)
<br>
NOTE 1: The introduction of the dominant reactive term involves replacing the matrix associated with the discretization of the bilinear form of the reactive term with the "mass lumping" matrix.
<br>
NOTE 2: Corrective terms are also added in the case of a dominant convective term.
<br>


Example of instability correction by using of "mass lumping":
<br>
![immagine](https://github.com/user-attachments/assets/f4b02890-cc5e-4e66-a57d-43bd54b081aa)


## Non-Stationary Problem Case
The code provides an example of applying the method of lines:
<br>
- The Finite Element Method (FEM) is used to approximate the solution of the PDE in space.
- From the system of ordinary differential equations obtained, a time solution is applied using the Implicit Euler method.
<br>
WARNING: For course-related reasons, the mesh creation part has been assigned to a Delaunay Triangulation package, which provides the data structures for the nodes and edges of the model. This package IS NOT INCLUDED IN THIS REPOSITORY.



