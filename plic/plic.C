/*---------------------------------------------------------------------------*\
    =========                 |
    \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
    \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
    \\/     M anipulation  |
    -------------------------------------------------------------------------------
    License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

    Author
    Ashwin Raghavan

    \*---------------------------------------------------------------------------*/

//May all your coding and simulations work perfectly to bring you Success in all your work.

#include "plic.H"
#include "vector.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "surfaceInterpolate.H"
#include "plane.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(plic, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::plic::intfc_cell_reconstruct
(
    const vector& intfc_normal,
    const scalar& alphaRef,
    const cell& curCell,
    const faceList& faces,
    const pointField& points,
    cell& ph1_cell,
    cell& ph0_cell,
    faceList& fcs,
    pointField& pts,
    face& intfc_face,
    labelList& ph1_fcLbls,
    labelList& ph0_fcLbls
)
{
    scalar nx = intfc_normal.x();
    scalar ny = intfc_normal.y();
    scalar nz = intfc_normal.z();

    if(debugIR2_)
    {
        Info<< "Interface normal: (" << nx << " " << ny << " " << nz << ")" << endl;
        Info<< "Determining initial guesses for plane constant..." << endl;
        Info<< "Candidate initial planes passing through cell points" << endl;
        Info<< "----------------------------------------------------------" << endl;
        Info<< "Point index            Point            Plane constant    " << endl;
        Info<< "----------------------------------------------------------" << endl;
    }

    const pointField& cellPts = curCell.points(faces, points);
    point curPt = cellPts[0];
    scalar a = -(nx*curPt.x() + ny*curPt.y() + nz*curPt.z());
    if(debugIR2_)
    {
        Info<< "    " << 0 << "          " << curPt << "         " << a << endl; 
    }
    scalar b = a;
    scalarList pln_coeffs(4);
    pln_coeffs[0] = nx; pln_coeffs[1] = ny; pln_coeffs[2] = nz;
    
    for(label pointI=1; pointI<cellPts.size(); pointI++)
    {
        curPt = cellPts[pointI];
        scalar curPt_pln_const = -(nx*curPt.x() + ny*curPt.y() + nz*curPt.z());
     
        if(debugIR2_)
        {
            Info<< "    " << pointI << "          " << curPt << "         " << curPt_pln_const << endl; 
        }
        
        if(curPt_pln_const < a)
        {            
            a = curPt_pln_const;
        }
        if(curPt_pln_const > b)
        {
            b = curPt_pln_const;
        }
    }

    if(debugIR2_)
    {
        Info<< "----------------------------------------------------------" << endl << endl;
    }
    
    scalar cellVol = curCell.mag(points, faces);

    if(debugIR2_)
    {
        Info<< "Cell volume  =   " << cellVol << endl << endl;
    }

    if(debugIR2_)
    {
        Info<< "------------ Initializing Brent iterations --------------" << endl;
    }

    if(debugIR2_)
    {
        Info<< "---------------------------------------------------------" << endl;
        Info<< "Clipping cell with plane with constant = " << a << endl << endl;
    }
    pln_coeffs[3] = a;
    Plane pln_a(pln_coeffs);    
    plane_cell_intersect(pln_a, curCell, faces, points, ph1_cell, ph0_cell, fcs, pts, intfc_face, ph1_fcLbls, ph0_fcLbls);    
    scalar fa = ph1_cell.mag(pts, fcs)/cellVol - alphaRef;
    if(debugIR2_)
    {        
        Info<< "Cell alpha1 error for plane  = " << fa << endl;
        Info<< "---------------------------------------------------------" << endl << endl;
    }
    
    if(debugIR2_)
    {
        Info<< "---------------------------------------------------------" << endl;
        Info<< "Clipping cell with plane with constant = " << b << endl << endl;
    }
    pln_coeffs[3] = b;
    Plane pln_b(pln_coeffs);
    subCells_reset(ph1_cell, ph0_cell, fcs, pts, intfc_face);
    plane_cell_intersect(pln_b, curCell, faces, points, ph1_cell, ph0_cell, fcs, pts, intfc_face, ph1_fcLbls, ph0_fcLbls);    
    scalar fb = ph1_cell.mag(pts, fcs)/cellVol - alphaRef;
    if(debugIR2_)
    {        
        Info<< "Cell alpha1 error for plane  = " << fb << endl;
        Info<< "---------------------------------------------------------" << endl << endl;
    }
    
    if(mag(fa) < mag(fb))
    {
        swap_scalar(a, b);
        swap_scalar(fa, fb);
    }

    scalar bOld = a; scalar fbOld = fa;
    scalar b_cpy = b; scalar fb_cpy = fb;
    label nIters = 0;
    scalar s1 = 0; scalar s2 = 0;

    if(debugIR2_)
    {
        Info<< "------------ Starting Brent iterations --------------" << endl << endl;
    }

    do
    {
        if(debugIR2_)
        {
            Info<< "-------------------------------------------------" << endl;
            Info<< "Brent iteration no:  " << nIters+1 << endl;
            Info<< "b = " << b << "    bOld = " << bOld << "    a = " << a << endl;
            Info<< "f(b) = " << fb << "    f(bOld) = " << fbOld << "    f(a) = " << fa << endl;
        }
        b_cpy = b; fb_cpy = fb;

        if(mag(fb - fbOld) < ALPHA_TOL)
        {
            b = 0.5*(b + a);
        }
        else
        {
            s1 = b - fb*(b - bOld)/(fb - fbOld);
            s2 = 0.5*(b + a);                    
            if(mag(b-s1) < mag(b-s2))
                b = s1;
            else
                b = s2;
        }
     
        if(debugIR2_)
        {
            Info<< "Secant soln. s1 = " << s1 << "   Bisection soln. s2 = " << s2 << endl; 
            Info<< "bNew  =  " << b << endl << endl;
            Info<< "---------------------------------------------------------" << endl;
            Info<< "Clipping cell with plane with constant = " << b << endl << endl;
        }
        pln_coeffs[3] = b;
        Plane pln(pln_coeffs);
        subCells_reset(ph1_cell, ph0_cell, fcs, pts, intfc_face);
        plane_cell_intersect(pln, curCell, faces, points, ph1_cell, ph0_cell, fcs, pts, intfc_face, ph1_fcLbls, ph0_fcLbls);    
        fb = ph1_cell.mag(pts, fcs)/cellVol - alphaRef;
        if(debugIR2_)
        {        
            Info<< "Cell alpha1 error for plane  = " << fb << endl;
            Info<< "---------------------------------------------------------" << endl << endl;
        }

        bOld = b_cpy; fbOld = fb_cpy;
        if(fb*fa > 0)
        {
            a = bOld;
            fa = fbOld;
        }
        if(debugIR2_)
        {
            Info<< "bNew = " << b << "    a = " << a << endl; 
            Info<< "f(bNew) = " << fb << "    f(a) = " << fa << endl; 
        }

        if(mag(fa) < mag(fb))
        {
            swap_scalar(a, b);
            swap_scalar(fa, fb);
        }
        if(debugIR2_)
        {
            Info<< "Swapping a and b if |f(a)| < |f(b)|" << endl;
            Info<< "bNew = " << b << "    a = " << a << endl; 
            Info<< "f(bNew) = " << fb << "    f(a) = " << fa << endl; 
            Info<< "-------------------------------------------------" << endl << endl;
        }
        
        nIters++;

    }while(mag(fb) > ALPHA_TOL && nIters < MAX_BRENT_ITERS);

    if(debugIR_)
    {
        Info<< "No. of Brent iterations:  " << nIters << endl << endl;
    }

    brent_iters_tmp_ = nIters;
    brent_err_tmp_ = mag(fb);
}
    

void Foam::plic::plane_cell_intersect
(
    const Plane& intfc_pln,
    const cell& curCell,
    const faceList& faces,
    const pointField& points,
    cell& ph1_cell,
    cell& ph0_cell,
    faceList& fcs,
    pointField& pts,
    face& intfc_face,
    labelList& ph1_fcLbls,
    labelList& ph0_fcLbls
)
{
    if(debugIR2_)
    {
        Info<< "==============================================================" << endl;
        Info<< "        Clipping Cell with given Plane        " << endl;
        Info<< "--------------------------------------------------------------" << endl << endl;
        Info<< "Current Cell:" << endl;
        Foam::plicFuncs::display_cell(curCell, faces, points);
        Info<< "Clipping Plane:" << endl;
        Foam::plicFuncs::display_Plane(intfc_pln);
    }

    labelList intfc_pts_lbls(0);

    bool intfc_on_cell_face = false;
    label intfc_cell_face_lbl;

    // Loop over all faces of the cell
    forAll(curCell, faceI)
    {                
        face curFace = faces[curCell[faceI]];
        labelList intfc_pts_lbls_curFace(0);

        if(debugIR2_)
        {
            Info<< "==============================================================" << endl;
            Info<< "        Clipping face " << faceI  << endl;
            Info<< "--------------------------------------------------------------" << endl << endl;
            Info<< "Current face:" << endl;
            Foam::plicFuncs::display_face(curFace, points);
        }

        // Distance of face points from Plane
        scalarField d_curPts(curFace.size());
        if(debugIR2_)
        {
            Info<< "--------------------------------------------------------------" << endl;
            Info<< "        Point                          Dist from plane        " << endl;
            Info<< "--------------------------------------------------------------" << endl;
        }
        forAll(curFace, pointI)
        {
            d_curPts[pointI] = calc_signed_dist(intfc_pln, points[curFace[pointI]]);
            if(debugIR2_)
            {
                Info<< "      " << points[curFace[pointI]] << "                    " << d_curPts[pointI] << endl;
            }
        }
        if(debugIR2_)
        {
            Info<< "--------------------------------------------------------------" << endl << endl;
        }        

        // First check if pure face
        bool pureFace = true;
        for(label i=1; i<curFace.size(); i++)
        {
            if(!(d_curPts[i]*d_curPts[0] > 0))
            {
                pureFace = false;
                break;
            }               
        }

        for(label i=0; i<curFace.size(); i++)
        {
            if(mag(d_curPts[i]) < DIST_TOL)
            {
                pureFace = false;
                break;
            }
        }

        label nIntfcPts = 0;
        for(label i=0; i<curFace.size(); i++)
        {
            if(mag(d_curPts[i]) < DIST_TOL)
            {
                nIntfcPts++;
            }
        }

        if(nIntfcPts >= 3)
        {
            intfc_on_cell_face = true;
            intfc_cell_face_lbl = faceI;
            for(label i=0; i<curFace.size(); i++)
            {
                d_curPts[i] = 0.1*DIST_TOL;            
            }
            pureFace = false;
        }

        if(debugIR2_)
        {
            Info<< "Pure phase indicator for face:  " << pureFace << endl << endl;
        }

        // If a pure face, assign the face to appropriate phase
        if(pureFace)
        { 
            face tFace(0);

            if(d_curPts[0] > 0)
            {
                for(label pointI=0; pointI<curFace.size(); pointI++)
                {
                    tFace.append(append_point_if_new(pts,points[curFace[pointI]]));                    
                }
                fcs.append(tFace);
                ph1_cell.append(fcs.size()-1);
                ph1_fcLbls[faceI] = fcs.size()-1;
                ph0_fcLbls[faceI] = -1;

                if(debugIR2_)
                {
                    Info<< "Face appended to phase 1 subcell: " << endl;
                    Foam::plicFuncs::display_face(tFace, pts);
                }
            }
            else if(d_curPts[0] < 0)
            {
                for(label pointI=0; pointI<curFace.size(); pointI++)
                {                    
                    tFace.append(append_point_if_new(pts,points[curFace[pointI]]));
                }

                tFace.collapse();
                fcs.append(tFace);
                ph0_cell.append(fcs.size()-1);
                ph0_fcLbls[faceI] = fcs.size()-1;
                ph1_fcLbls[faceI] = -1;

                if(debugIR2_)
                {
                    Info<< "Face appended to phase 0 subcell: " << endl;
                    Foam::plicFuncs::display_face(tFace, pts);
                }
            }            
        } // end of loop for pure face
        else // i.e. if mixed face
        {
            face tFace_ph1(0); 
            face tFace_ph0(0);
            
            // loop over all edges of face
            for(label edgeI=0; edgeI<curFace.size(); edgeI++)
            {
                if(debugIR2_)
                {
                    Info<< "---------------------------------------------------------------------" << endl;
                    Info<< "Edge " << edgeI << " intersection with Plane" << endl;
                    Info<< "---------------------------------------------------------------------" << endl; 
                }
                label startPt_lbl = edgeI;
                label endPt_lbl;
                if(edgeI < curFace.size()-1)
                {
                    endPt_lbl = edgeI + 1;
                }
                else
                {
                    endPt_lbl = 0;
                }
                
                point startPt = points[curFace[startPt_lbl]];
                point endPt = points[curFace[endPt_lbl]];
                scalar d1 = d_curPts[startPt_lbl];
                scalar d2 = d_curPts[endPt_lbl];

                if(debugIR2_)
                {
                    Info<< "start pt:  " << startPt << "      start pt dist:  " << d1 << endl;
                    Info<< "end pt:  " << endPt << "      end pt dist:  " << d2 << endl;
                }

                //start ifs
                if(mag(d1) < DIST_TOL || mag(d2) < DIST_TOL)
                {
                    if(mag(d1) < DIST_TOL && mag(d2) < DIST_TOL)
                    {
                        if(debugIR2_)
                        {
                            Info<< "d1=0 & d2=0; both pts=>ph1; both pts=>ph0; both pts=>intfc face pts" << endl;
                        }
                        tFace_ph1.append(append_point_if_new(pts, startPt));
                        intfc_pts_lbls_curFace.append(tFace_ph1.last());
                        tFace_ph1.append(append_point_if_new(pts, endPt));
                        intfc_pts_lbls_curFace.append(tFace_ph1.last());
                        tFace_ph0.append(append_point_if_new(pts, startPt));
                        tFace_ph0.append(append_point_if_new(pts, endPt));
                    }
                    else if(mag(d1) < DIST_TOL)
                    {
                        if(d2 > 0)
                        {
                            if(debugIR2_)
                            {
                                Info<< "d1=0 & d2>0; both pts=>ph1; start pt=>ph0; start pt=>intfc face pts" << endl;
                            }
                            tFace_ph0.append(append_point_if_new(pts, startPt));             
                            tFace_ph1.append(append_point_if_new(pts, startPt));
                            intfc_pts_lbls_curFace.append(tFace_ph1.last());
                            tFace_ph1.append(append_point_if_new(pts, endPt));
                        }
                        else
                        {
                            if(debugIR2_)
                            {
                                Info<< "d1=0 & d2<0; start pt=>ph1; both pts=>ph0; start pt=>intfc face pts" << endl;
                            }
                            tFace_ph1.append(append_point_if_new(pts, startPt));             
                            intfc_pts_lbls_curFace.append(tFace_ph1.last());
                            tFace_ph0.append(append_point_if_new(pts, startPt));
                            tFace_ph0.append(append_point_if_new(pts, endPt));
                        }
                    }
                    else
                    {
                        if(d1 > 0)
                        {
                            if(debugIR2_)
                            {
                                Info<< "d2=0 & d1>0; both pts=>ph1; end pt=>ph0; end pt=>intfc face pts" << endl;
                            }
                            tFace_ph1.append(append_point_if_new(pts, startPt));             
                            tFace_ph1.append(append_point_if_new(pts, endPt));
                            intfc_pts_lbls_curFace.append(tFace_ph1.last());
                            tFace_ph0.append(append_point_if_new(pts, endPt));
                        }
                        else
                        {                       
                            if(debugIR2_)
                            {
                                Info<< "d2=0 & d1<0; both pts=>ph0; end pt=>ph1; end pt=>intfc face pts" << endl;
                            }
                            tFace_ph0.append(append_point_if_new(pts, startPt));
                            tFace_ph0.append(append_point_if_new(pts, endPt));
                            tFace_ph1.append(append_point_if_new(pts, endPt));               
                            intfc_pts_lbls_curFace.append(tFace_ph1.last());
                        }
                    }
                }
                else
                {
                    if(d1*d2 > 0)
                    {                    
                        if(d1 > 0)
                        {                
                            if(debugIR2_)
                            {
                                Info<< "d1*d2>0 & d1>0; both pts=>ph1" << endl;
                            }
                            tFace_ph1.append(append_point_if_new(pts, startPt));
                            tFace_ph1.append(append_point_if_new(pts, endPt));
                        }
                        else
                        {   
                            if(debugIR2_)
                            {
                                Info<< "d1*d2>0 & d1<0; both pts=>ph0" << endl;
                            }
                            tFace_ph0.append(append_point_if_new(pts, startPt));
                            tFace_ph0.append(append_point_if_new(pts, endPt));
                        }
                    }
                    else if(d1*d2 < 0)
                    {
                        point intfc_pt = plane_edge_intersect(startPt, endPt, d1, d2);

                        if(debugIR2_)
                        {
                            Info<< "Interface point:  " << intfc_pt << endl;
                        }

                        if(d1 > 0)
                        {
                            if(debugIR2_)
                            {
                                Info<< "d1*d2<0 & d1>0; start pt,intfc pt=>ph1; intfc pt,end pt=>ph0; intfc pt=>intfc face pts" << endl;
                            }
                            tFace_ph1.append(append_point_if_new(pts, startPt));
                            tFace_ph1.append(append_point_if_new(pts, intfc_pt));
                            intfc_pts_lbls_curFace.append(tFace_ph1.last());
                            tFace_ph0.append(append_point_if_new(pts, intfc_pt));
                            tFace_ph0.append(append_point_if_new(pts, endPt));
                        }
                        else
                        {
                            if(debugIR2_)
                            {
                                Info<< "d1*d2<0 & d1<0; start pt,intfc pt=>ph0; intfc pt,end pt=>ph1; intfc pt=>intfc face pts" << endl;
                            }
                            tFace_ph0.append(append_point_if_new(pts, startPt));
                            tFace_ph0.append(append_point_if_new(pts, intfc_pt));     
                            tFace_ph1.append(append_point_if_new(pts, intfc_pt));
                            intfc_pts_lbls_curFace.append(tFace_ph1.last());
                            tFace_ph1.append(append_point_if_new(pts, endPt));
                        }
                    }
                }                
                //end ifs

                if(debugIR2_)
                {
                    Info<< "-----------------------------------------------------------------------" << endl;
                    Info<< "Done edge " << edgeI << " intersection with Plane" << endl;
                    Info<< "-----------------------------------------------------------------------" << endl << endl; 
                }
            } // end of loop over edges

            tFace_ph1.collapse();
            tFace_ph0.collapse();

            if(debugIR2_)
            {
                Info<< "Interface points list for current face: " << endl;
                Foam::plicFuncs::display_labelList(intfc_pts_lbls_curFace);
            }

            remove_duplicateLabels(intfc_pts_lbls_curFace);

            if(debugIR2_)
            {
                Info<< "Interface points list for current face without duplicates: " << endl;
                Foam::plicFuncs::display_labelList(intfc_pts_lbls_curFace);
            }
            
            bool pairRepeated = false;
            if(intfc_pts_lbls_curFace.size() == 2)
            {
                for(label lblI=0; lblI<intfc_pts_lbls.size(); lblI++)
                {
                    if(intfc_pts_lbls_curFace[0] == intfc_pts_lbls[lblI])
                    {
                        if(lblI % 2)
                        {
                            if(intfc_pts_lbls_curFace[1] == intfc_pts_lbls[lblI-1])
                            {
                                pairRepeated = true;
                                break;
                            }
                        }
                        else
                        {
                            if(intfc_pts_lbls_curFace[1] == intfc_pts_lbls[lblI+1])
                            {
                                pairRepeated = true;
                                break;
                            }
                        }
                    }
                }
                
                if(!pairRepeated)
                {                
                    intfc_pts_lbls.append(intfc_pts_lbls_curFace[0]);
                    intfc_pts_lbls.append(intfc_pts_lbls_curFace[1]);
                }
            }

            if(debugIR2_)
            {
                Info<< "Interface points pair repeated?:  " << pairRepeated << nl << nl
                    << "Point field: " << endl;
                Foam::plicFuncs::display_pointField(pts);
                Info<< "Phase 1 sub-face: " << endl;
                Foam::plicFuncs::display_face(tFace_ph1, pts);
                Info<< "Phase 0 sub-face: " << endl;
                Foam::plicFuncs::display_face(tFace_ph0, pts);
                Info<< "Interface points list: " << endl;
                Foam::plicFuncs::display_labelList(intfc_pts_lbls);
            }

            if(tFace_ph1 == tFace_ph0)
            {
                ph1_fcLbls[faceI] = -2;
                ph0_fcLbls[faceI] = -2;
            }
            else
            {
                if(tFace_ph1.size() > 2)
                {
                    fcs.append(tFace_ph1);
                    ph1_cell.append(fcs.size()-1);
                    ph1_fcLbls[faceI] = fcs.size()-1;
                }
                else
                {
                    ph1_fcLbls[faceI] = -1;
                }
                if(tFace_ph0.size() > 2)
                {                                
                    fcs.append(tFace_ph0);
                    ph0_cell.append(fcs.size()-1);
                    ph0_fcLbls[faceI] = fcs.size()-1;                
                }
                else
                {
                    ph0_fcLbls[faceI] = -1;
                }
            }
        } // end of mixed face
        if(debugIR2_)
        {
            Info<< "--------------------------------------------------------------" << endl;            
            Info<< "        Done Clipping face " << faceI  << endl;
            Info<< "==============================================================" << endl << endl;
        }   
    } // end of loop over all cell faces        

    bool intfc_empty = true;
    if(intfc_pts_lbls.size())
    {
        intfc_empty = false;
    }

    if(intfc_on_cell_face)
    {
        if(debugIR2_)
        {
            Info<< "Cell face " << intfc_cell_face_lbl << " is the interface face..." << nl
                << "Skipping interface face construction..." << nl << endl;
        }

        face intfc_cell_face = faces[curCell[intfc_cell_face_lbl]];
        for(label pointI=0; pointI<intfc_cell_face.size(); pointI++)
        {            
            intfc_face.append(append_point_if_new(pts, points[intfc_cell_face[pointI]]));
        }
    }

    if(!intfc_empty && !intfc_on_cell_face)
    {
        if(debugIR2_)
        {
            Info<< "==============================================================" << endl;
            Info<< "                  Constructing interface face                 " << endl;
            Info<< "--------------------------------------------------------------" << endl << endl;
            Info<< "Interface point labels list: " << endl;
            Foam::plicFuncs::display_labelList(intfc_pts_lbls);
        }

        // reconstructing the cap face i.e. the interface face
        if(debugIR2_)
        {
            Info<< "--------------------------------------------------------------" << endl;
            Info<< "                      Iteration: 0" << endl;
            Info<< "--------------------------------------------------------------" << endl;
        }

        intfc_face.append(intfc_pts_lbls[0]);
        intfc_face.append(intfc_pts_lbls[1]);
        remove_elems(intfc_pts_lbls, 0, 1);

        if(debugIR2_)
        {
            Info<< "Interface point labels list: " << endl;
            Foam::plicFuncs::display_labelList(intfc_pts_lbls);
            Info<< "Interface face: " << endl;
            Foam::plicFuncs::display_face(intfc_face, pts);
            Info<< "--------------------------------------------------------------" << endl;
            Info<< "                    Done Iteration: 0" << endl;
            Info<< "--------------------------------------------------------------" << endl << endl;
        }

        label nIters = 1;
        do
        {
            if(debugIR2_)
            {
                Info<< "--------------------------------------------------------------" << endl;
                Info<< "                      Iteration: " << nIters << endl;
                Info<< "--------------------------------------------------------------" << endl;
            }
            for(int i=0; i < intfc_pts_lbls.size(); i++)
            {
                if(intfc_pts_lbls[i] == intfc_face.last())
                {
                    if(debugIR2_)
                    {
                        Info<< "Matching pt lbl found at i = " << i << endl;
                    }
                    if(i % 2)
                    {
                        if(intfc_pts_lbls[i-1] != intfc_face[intfc_face.size()-2])
                        {
                            intfc_face.append(intfc_pts_lbls[i-1]);
                            if(debugIR2_)
                            {
                                Info<< "i is odd & i-1 lbl != intfc_face(last-1); remove i-1,i lbls frm intfc pts lbls; i-1 lbl=>intfc_face" << endl;
                            }
                        }
                        remove_elems(intfc_pts_lbls, i-1, i);
                        if(debugIR2_)
                        {
                            Info<< "i is odd & i-1 lbl = intfc_face(last-1); remove i-1,i lbls frm intfc pts lbls" << endl;
                        }
                    }
                    else
                    {
                        if(intfc_pts_lbls[i+1] != intfc_face[intfc_face.size()-2])
                        {
                            intfc_face.append(intfc_pts_lbls[i+1]);
                            if(debugIR2_)
                            {
                                Info<< "i is even & i+1 lbl != intfc_face(last-1); remove i,i+1 lbls frm intfc pts lbls; i+1 lbl=>intfc_face" << endl;
                            }
                        }
                        remove_elems(intfc_pts_lbls, i, i+1);
                        if(debugIR2_)
                        {
                            Info<< "i is even & i+1 lbl = intfc_face(last-1); remove i,i+1 lbls frm intfc pts lbls" << endl;
                        }
                    }
                    break;
                }
            }
            if(debugIR2_)
            {
                Info<< "Interface point labels list: " << endl;
                Foam::plicFuncs::display_labelList(intfc_pts_lbls);
                Info<< "Interface face: " << endl;
                Foam::plicFuncs::display_face(intfc_face, pts);
                Info<< "--------------------------------------------------------------" << endl;
                Info<< "                    Done Iteration: " << nIters << endl;
                Info<< "--------------------------------------------------------------" << endl << endl;
            }
            nIters++;

            if(nIters > 20)
            {
                FatalErrorIn("void plic::plane_cell_intersect(const Plane&,const cell&,const faceList&,const pointField&,cell&,cell&,faceList&,pointField&,face&,labelList&,labelList&)")
                    << "Intfc face construction failed"
                        << abort(FatalError);
            }

        }while(intfc_pts_lbls.size() > 0);

        if(debugIR2_)
        {
            Info<< "Constructed interface face: " << endl;
            Foam::plicFuncs::display_face(intfc_face, pts);
        }

        intfc_face.collapse();

        if(debugIR2_)
        {
            Info<< "Collapsed interface face: " << endl;
            Foam::plicFuncs::display_face(intfc_face, pts);
            Info<< "Interface face normal: " << intfc_face.normal(pts) << endl;
        }
        
        if(intfc_face.size() > 2)
        {
            if((intfc_face.normal(pts) & intfc_pln.normal()) < 0)
            {
                intfc_face.flip();
                if(debugIR2_)
                {
                    Info<< "Reversed interface face orientation ..." << endl << endl;
                }
            }
        }

        if(debugIR2_)
        {
            Info<< "--------------------------------------------------------------" << endl;
            Info<< "  Done constructing interface face      No. of iterations: " << nIters << endl;
            Info<< "==============================================================" << endl << endl;
        }
        if(intfc_face.size() > 2)
        {
            fcs.append(intfc_face);
            ph1_cell.append(fcs.size()-1);
            ph0_cell.append(fcs.size()-1);        
        }// add intfc_face to both subcells only if it has more than 2 pts
    } // interface face is constructed only if intfc_pts_lbls was not empty
    
    if(debugIR2_)
    {
        Info<< "Interface face: " << endl;
        Foam::plicFuncs::display_face(intfc_face, pts);
        Info<< "Phase 1 sub-cell: " << endl;
        Foam::plicFuncs::display_cell(ph1_cell, fcs, pts);
        Info<< "Phase 0 sub-cell: " << endl;
        Foam::plicFuncs::display_cell(ph0_cell, fcs, pts);
        Info<< "--------------------------------------------------------------" << endl;    
        Info<< "        Done Clipping Cell with given Plane        " << endl;        
        Info<< "==============================================================" << endl << endl;
    }
}


Foam::scalar Foam::plic::calc_signed_dist
(
    const Plane& pln,
    const point& p
)
{
    return ((p - pln.refPoint()) & pln.normal());
}


Foam::point Foam::plic::plane_edge_intersect
(
    const point& x1,
    const point& x2,
    const scalar& d1,
    const scalar& d2
)
{
    scalar lambda = d1/(d2 - d1);
    return (x1 - lambda*(x2 - x1));
}


Foam::label Foam::plic::append_point_if_new
(
    pointField& Pts,
    const point& pt
)
{
    if(debugIR2_)
    {
        Info<< "---------------------------------------------------------------" << endl;
        Info<< "Checking if new point..." << endl;
    }

    //bool newPt = true;
    for(label i=0; i<Pts.size(); i++)
    {
        if(debugIR2_)
        {
            Info<< "Candidate point: " << pt << "    pts[" << i << "]: " << Pts[i] << endl;
            Info<< "Candidate point == pts[i]?: " << isEq_pt(pt, Pts[i]) << endl;
        }

        if(isEq_pt(pt, Pts[i]))
        {
            //newPt = false;
            if(debugIR2_)
            {        
                Info<< "Not a new point..." << endl;
                Info<< "---------------------------------------------------------------" << endl << endl;
            }
            return i;
        }
    }
    if(debugIR2_)
    {        
        Info<< "New point appended to pointField" << endl;
        Info<< "---------------------------------------------------------------" << endl << endl;
    }
    //if(newPt)
    //{
    Pts.append(pt);
    return (Pts.size() - 1);
    //}
}


void Foam::plic::remove_elems
(
    labelList& lbls,
    const label a,
    const label b
)
{
    labelList lbls_copy(lbls);
    lbls.resize(lbls.size()-2);
    label last = 0;
    for(label i=0; i<lbls_copy.size(); i++)
    {
        if(i!=a && i!=b)
        {
            lbls[last] = lbls_copy[i];
            last++;
        }
    }
}


void Foam::plic::remove_duplicateLabels
(
    labelList& lbls
)
{
    labelList lbls_copy(lbls);
    lbls.resize(1);
    lbls[0] = lbls_copy[0];
    bool duplicate = false;
    for(label i=1; i<lbls_copy.size(); i++)
    {
        duplicate = false;
        for(label j=0; j<lbls.size(); j++)
        {
            if(lbls_copy[i] == lbls[j])
            {
                duplicate = true;
                break;
            }
        }
        if(!duplicate)
        {
            lbls.append(lbls_copy[i]);
        }
    }
}


void Foam::plic::subCells_reset
(
    cell& ph1_cell, 
    cell& ph0_cell, 
    faceList& fcs, 
    pointField& pts, 
    face& intfc_face
)
{
    ph1_cell.setSize(0);
    ph0_cell.setSize(0);
    fcs.setSize(0);
    pts.setSize(0);
    intfc_face.setSize(0);        
}


void Foam::plic::swap_scalar
(
    scalar& a,
    scalar& b
)
{
    scalar tmp = a;
    a = b;
    b = tmp;
}


void Foam::plic::calcFaceFluxTets
(
    const label curFcI,
    const scalar& curPhi,
    List<tetPoints>& tets,
    List<scalar>& tetVols
)
{
    const faceList& faces = mesh_.faces();
    const pointField& points = mesh_.points();
    const labelListList& pointCells = mesh_.pointCells();
    const vectorField& UCells = U_.internalField();
    const labelList& own = mesh_.owner();
    const labelList& nei = mesh_.neighbour();

    face curFc = faces[curFcI];

    // Construct the face flux polyhedron by Lagrangian backtracking of 
    // face vertices. Calculate the centroid of backtracked face so
    // as to obtain a face flux polyhedron consistent with mass
    // conservation. Finally, decompose the polyhedron into 20
    // tetrahedrons which will give the correct sign of their 
    // contribution to the face flux based on their placement w.r.t
    // the face, that is, inflow or outflow tet

    pointField pts(8);
    faceList fcs(6);
    pointField fcCtrs(6);    

    scalar dt = mesh().time().deltaT().value();

    if(debugF_)
    {
        Info<< "Current time step:  " << dt << nl << nl
            << "--------------------------------------------------------------------" << nl
            << "            Lagrangian backtracking of face vertices" << nl
            << "--------------------------------------------------------------------" << nl << endl;
    }

    vector U_fcCtr = vector::zero;
    label nFcPts = 0;

    for(label ptI=0; ptI<curFc.size(); ptI++)
    {
        label curPtI = curFc[ptI];
        point curPt = points[curPtI];
        pts[ptI] = curPt;
        pts[ptI+4] = curPt - U_pts_[curPtI]*dt;

        U_fcCtr += U_pts_[curPtI];
        nFcPts++;        
        
        if(debugF_)
        {
            Info<< "vertex " << ptI << ":  " << pts[ptI] << nl
                << "velocity vector:  " << U_pts_[curPtI] << nl
                << "backtrack vector:  " << U_pts_[curPtI]*dt << nl
                << "backtrack vertex " << (ptI+4) << ":  " << pts[ptI+4] << nl << endl;
        }

        const labelList& pcp = pointCells[curPtI];

        if(debugF_)
        {
            forAll(pcp, cI)
            {
                if(pcp[cI]<mesh_.nCells())
                {
                    Info<< "cell " << pcp[cI] << " vel = " << UCells[pcp[cI]] << endl;
                }
            }
            Info<< endl;
        }
    }
    U_fcCtr /= nFcPts;
    vector curFcArea = Sf_[curFcI];
    scalar curFcCtrPhi = U_fcCtr & curFcArea;   

    if(debugF_)
    {
        Info<< "Face flux: " << curPhi << nl
            << "Face vel. based on flux: " << curPhi*Sf_[curFcI]/magSf_[curFcI]/magSf_[curFcI] << nl
            << "Face vel. based on cells: " << (UCells[own[curFcI]]+UCells[nei[curFcI]])/2.0 << nl
            << "Face vel. based on points: " << U_fcCtr << nl
            << "Face own vel.: " << UCells[own[curFcI]] << "  Face nei vel.: " << UCells[nei[curFcI]] << nl
            << "Flux based on face ctr vel.: " << curFcCtrPhi << nl
            << "Flux polyhedron vertex points:" << endl;
        Foam::plicFuncs::display_pointField(pts);
        Info<< "--------------------------------------------------------------------" << nl
            << "         Done Lagrangian backtracking of face vertices" << nl
            << "--------------------------------------------------------------------" << nl << endl;
    }

    fcs[0] = face_from_lbls(0,1,2,3);
    fcs[1] = face_from_lbls(0,4,5,1);
    fcs[2] = face_from_lbls(1,5,6,2);
    fcs[3] = face_from_lbls(2,6,7,3);
    fcs[4] = face_from_lbls(3,7,4,0);
    fcs[5] = face_from_lbls(4,7,6,5);

    for(label fcI=0; fcI<fcCtrs.size(); fcI++)
    {
        fcCtrs[fcI] = fcs[fcI].centre(pts);

        if(debugF_)
        {
            Info<< "Face " << fcI << " of flux polyhedron: " << endl;
            Foam::plicFuncs::display_face(fcs[fcI], pts);
            Info<< "Face centroid:  " << fcCtrs[fcI] << nl << endl;
        }
    }

    if(debugF_)
    {
        Info<< nl
            << "Constructing first 16 flux tetrahedra..." << nl << endl;
    }
    
    point a = fcCtrs[0];
    scalar V16 = 0;
    label tetI = 0;

    for(label fcI=1; fcI<(fcs.size()-1); fcI++)
    {
        face curFace = fcs[fcI];
        label nPts = curFace.size();
        point b = fcCtrs[fcI];
        for(label ptI=0; ptI<nPts; ptI++)
        {
            point c = pts[curFace[ptI]];
            point d = pts[curFace[(ptI + 1) % nPts]];
            tets[tetI] = tetPoints(a, b, c, d);
            tetVols[tetI] = tets[tetI].tet().mag();

            if(debugF_)
            {
                Info<< "Tetrahedron " << tetI << nl
                    << "Points: " << a << "  " << b << "  " << c << "  " << d << nl
                    << "Volume: " << tetVols[tetI] << nl << endl;
            }

            V16 += tetVols[tetI];
            tetI++;
        } 
    }    

    scalar L02 = fcs[0].mag(pts);
    //scalar L02 = L0*L0;
    //scalar L03 = L0*L0*L0;
    scalar Vf_curFc = curPhi*dt/L02;
    V16 /= L02;
    //vector niHat = (fcCtrs[5] - fcCtrs[0]);
    face curFace = fcs[5];
    vector niHat = (pts[curFace[0]] - fcCtrs[5]) ^ (pts[curFace[1]] - fcCtrs[5]);
    niHat /= (mag(niHat) + VSMALL);

    if(debugF_)
    {
        Info<< "Calculating face centre of last face..." << nl
            << "Current face flux: " << curPhi << nl
            << "Normalization factor L02: " << L02 << nl
            << "Total flux volume: " << Vf_curFc*L02 << nl
            << "Vol of 1st 16 tets: " << V16*L02 << nl
            << "Vector from 1st face to last face centre: " << niHat << endl;
    }

    vector B = vector::zero;
    //face curFace = fcs[5];
    label nPts = curFace.size();
    for(label ptI=0; ptI<nPts; ptI++)
    {
        point c = pts[curFace[ptI]];
        point d = pts[curFace[(ptI + 1) % nPts]];
        B += ((c - a) ^ (d - a))/L02;
    } 

    scalar den = (niHat & B);
    /*if(mag(den)<VSMALL)
        {
        den += VSMALL;
        }*/

    scalar lambda = 6.0*(Vf_curFc - V16)/(mag(den) + VSMALL)*sign(den);
    point n = fcCtrs[0] + lambda*niHat;
    fcCtrs[5] = n;

    if(debugF_)
    {
        Info<< "Vector B: " << B*L02 << nl
            << "denominator: " << den*L02 << nl
            << "lambda: " << lambda << nl
            << "Last face centre calculated: " << n << nl << nl
            << "Constructing last 4 flux tetrahedra..." << nl << endl;
    }

    V16 *= L02;
    Vf_curFc *= L02;

    scalar V20 = V16;
    point b = fcCtrs[5];
    for(label ptI=0; ptI<nPts; ptI++)
    {
        point c = pts[curFace[ptI]];
        point d = pts[curFace[(ptI + 1) % nPts]];
        tets[tetI] = tetPoints(a, b, c, d);
        tetVols[tetI] = tets[tetI].tet().mag();
        V20 += tetVols[tetI];
        if(debugF_)
        {
            Info<< "Tetrahedron " << tetI << nl
                << "Points: " << a << "  " << b << "  " << c << "  " << d << nl
                << "Volume: " << tetVols[tetI] << nl << endl;
        }
        tetI++;
    }

    scalar Vf_err = mag(Vf_curFc - V20);

    if(debugF_)
    {
        Info<< "Face flux volume: " << Vf_curFc << " Total volume of face flux tets: " << V20 << nl
            << "Face flux tetrahedrons flux error: " << Vf_err << nl << endl;
    }

    // end of tetrahedral decomposition of face flux polygon volume
}


void Foam::plic::calcFaceFluxTets
(
    const label curFcI,
    const scalar& curPhi,
    List<tetPoints>& tets,
    List<scalar>& tetVols,
    const scalar& deltaT
)
{
    const faceList& faces = mesh_.faces();
    const pointField& points = mesh_.points();
    const labelListList& pointCells = mesh_.pointCells();
    const vectorField& UCells = U_.internalField();
    const labelList& own = mesh_.owner();
    const labelList& nei = mesh_.neighbour();

    face curFc = faces[curFcI];

    // Construct the face flux polyhedron by Lagrangian backtracking of 
    // face vertices. Calculate the centroid of backtracked face so
    // as to obtain a face flux polyhedron consistent with mass
    // conservation. Finally, decompose the polyhedron into 20
    // tetrahedrons which will give the correct sign of their 
    // contribution to the face flux based on their placement w.r.t
    // the face, that is, inflow or outflow tet

    pointField pts(8);
    faceList fcs(6);
    pointField fcCtrs(6);    

    scalar dt = deltaT;

    if(debugF_)
    {
        Info<< "Current time step:  " << dt << nl << nl
            << "--------------------------------------------------------------------" << nl
            << "            Lagrangian backtracking of face vertices" << nl
            << "--------------------------------------------------------------------" << nl << endl;
    }

    vector U_fcCtr = vector::zero;
    label nFcPts = 0;

    for(label ptI=0; ptI<curFc.size(); ptI++)
    {
        label curPtI = curFc[ptI];
        point curPt = points[curPtI];
        pts[ptI] = curPt;
        pts[ptI+4] = curPt - U_pts_[curPtI]*dt;

        U_fcCtr += U_pts_[curPtI];
        nFcPts++;        
        
        if(debugF_)
        {
            Info<< "vertex " << ptI << ":  " << pts[ptI] << nl
                << "velocity vector:  " << U_pts_[curPtI] << nl
                << "backtrack vector:  " << U_pts_[curPtI]*dt << nl
                << "backtrack vertex " << (ptI+4) << ":  " << pts[ptI+4] << nl << endl;
        }

        const labelList& pcp = pointCells[curPtI];

        if(debugF_)
        {
            forAll(pcp, cI)
            {
                if(pcp[cI]<mesh_.nCells())
                {
                    Info<< "cell " << pcp[cI] << " vel = " << UCells[pcp[cI]] << endl;
                }
            }
            Info<< endl;
        }
    }
    U_fcCtr /= nFcPts;
    vector curFcArea = Sf_[curFcI];
    scalar curFcCtrPhi = U_fcCtr & curFcArea;   

    if(debugF_)
    {
        Info<< "Face flux: " << curPhi << nl
            << "Face vel. based on flux: " << curPhi*Sf_[curFcI]/magSf_[curFcI]/magSf_[curFcI] << nl
            << "Face vel. based on cells: " << (UCells[own[curFcI]]+UCells[nei[curFcI]])/2.0 << nl
            << "Face vel. based on points: " << U_fcCtr << nl
            << "Face own vel.: " << UCells[own[curFcI]] << "  Face nei vel.: " << UCells[nei[curFcI]] << nl
            << "Flux based on face ctr vel.: " << curFcCtrPhi << nl
            << "Flux polyhedron vertex points:" << endl;
        Foam::plicFuncs::display_pointField(pts);
        Info<< "--------------------------------------------------------------------" << nl
            << "         Done Lagrangian backtracking of face vertices" << nl
            << "--------------------------------------------------------------------" << nl << endl;
    }

    fcs[0] = face_from_lbls(0,1,2,3);
    fcs[1] = face_from_lbls(0,4,5,1);
    fcs[2] = face_from_lbls(1,5,6,2);
    fcs[3] = face_from_lbls(2,6,7,3);
    fcs[4] = face_from_lbls(3,7,4,0);
    fcs[5] = face_from_lbls(4,7,6,5);

    for(label fcI=0; fcI<fcCtrs.size(); fcI++)
    {
        fcCtrs[fcI] = fcs[fcI].centre(pts);

        if(debugF_)
        {
            Info<< "Face " << fcI << " of flux polyhedron: " << endl;
            Foam::plicFuncs::display_face(fcs[fcI], pts);
            Info<< "Face centroid:  " << fcCtrs[fcI] << nl << endl;
        }
    }

    if(debugF_)
    {
        Info<< nl
            << "Constructing first 16 flux tetrahedra..." << nl << endl;
    }
    
    point a = fcCtrs[0];
    scalar V16 = 0;
    label tetI = 0;

    for(label fcI=1; fcI<(fcs.size()-1); fcI++)
    {
        face curFace = fcs[fcI];
        label nPts = curFace.size();
        point b = fcCtrs[fcI];
        for(label ptI=0; ptI<nPts; ptI++)
        {
            point c = pts[curFace[ptI]];
            point d = pts[curFace[(ptI + 1) % nPts]];
            tets[tetI] = tetPoints(a, b, c, d);
            tetVols[tetI] = tets[tetI].tet().mag();

            if(debugF_)
            {
                Info<< "Tetrahedron " << tetI << nl
                    << "Points: " << a << "  " << b << "  " << c << "  " << d << nl
                    << "Volume: " << tetVols[tetI] << nl << endl;
            }

            V16 += tetVols[tetI];
            tetI++;
        } 
    }    

    scalar L02 = fcs[0].mag(pts);
    //scalar L02 = L0*L0;
    //scalar L03 = L0*L0*L0;
    scalar Vf_curFc = curPhi*dt/L02;
    V16 /= L02;
    //vector niHat = (fcCtrs[5] - fcCtrs[0]);
    face curFace = fcs[5];
    vector niHat = (pts[curFace[0]] - fcCtrs[5]) ^ (pts[curFace[1]] - fcCtrs[5]);
    niHat /= (mag(niHat) + VSMALL);

    if(debugF_)
    {
        Info<< "Calculating face centre of last face..." << nl
            << "Current face flux: " << curPhi << nl
            << "Normalization factor L02: " << L02 << nl
            << "Total flux volume: " << Vf_curFc*L02 << nl
            << "Vol of 1st 16 tets: " << V16*L02 << nl
            << "Vector from 1st face to last face centre: " << niHat << endl;
    }

    vector B = vector::zero;
    //face curFace = fcs[5];
    label nPts = curFace.size();
    for(label ptI=0; ptI<nPts; ptI++)
    {
        point c = pts[curFace[ptI]];
        point d = pts[curFace[(ptI + 1) % nPts]];
        B += ((c - a) ^ (d - a))/L02;
    } 

    scalar den = (niHat & B);
    /*if(mag(den)<VSMALL)
        {
        den += VSMALL;
        }*/

    scalar lambda = 6.0*(Vf_curFc - V16)/(mag(den) + VSMALL)*sign(den);
    point n = fcCtrs[0] + lambda*niHat;
    fcCtrs[5] = n;

    if(debugF_)
    {
        Info<< "Vector B: " << B*L02 << nl
            << "denominator: " << den*L02 << nl
            << "lambda: " << lambda << nl
            << "Last face centre calculated: " << n << nl << nl
            << "Constructing last 4 flux tetrahedra..." << nl << endl;
    }

    V16 *= L02;
    Vf_curFc *= L02;

    scalar V20 = V16;
    point b = fcCtrs[5];
    for(label ptI=0; ptI<nPts; ptI++)
    {
        point c = pts[curFace[ptI]];
        point d = pts[curFace[(ptI + 1) % nPts]];
        tets[tetI] = tetPoints(a, b, c, d);
        tetVols[tetI] = tets[tetI].tet().mag();
        V20 += tetVols[tetI];
        if(debugF_)
        {
            Info<< "Tetrahedron " << tetI << nl
                << "Points: " << a << "  " << b << "  " << c << "  " << d << nl
                << "Volume: " << tetVols[tetI] << nl << endl;
        }
        tetI++;
    }

    scalar Vf_err = mag(Vf_curFc - V20);

    if(debugF_)
    {
        Info<< "Face flux volume: " << Vf_curFc << " Total volume of face flux tets: " << V20 << nl
            << "Face flux tetrahedrons flux error: " << Vf_err << nl << endl;
    }

    // end of tetrahedral decomposition of face flux polygon volume
}

/*
void Foam::plic::makeCorrVecFlatFld()
{
    label cellI, nCompact, nCells, nFlatFld;

    const fvMesh& mesh = mesh();
    const mapDistribute& map = faceStencil().map();

    nCells = mesh.nCells();
    
    nFlatFld = map.constructSize();
    corrVec_flatFld_.resize(nFlatFld);
    

    for(cellI=0; cellI<nCells; cellI++)
    {
        corrVec_flatFld_[cellI] = vector::zero;
    }

    forAll(alpha_ph1_.boundaryField(), patchI)
    {
        const polyPatch& pp = mesh.boundary()[patchI].patch();
        const fvPatchScalarField& pAlpha1 = alpha_ph1_.boundaryField()[patchI];

        nCompact = pp.start() - mesh.nInternalFaces() + nCells;

        if(isA<processorCyclicPolyPatch>(pp))
        {
            const vectorField& sep = pp.separation();

            forAll(pAlpha1, fcI)
            {
                corrVec_flatFld_[nCompact++] = sep[fcI];
            }
        }
        else
        {
            forAll(pAlpha1, fcI)
            {
                corrVec_flatFld_[nCompact++] = vector::zero;
            }
        }
    }

    map.distribute(corrVec_flatFld_);
}
    */

Foam::face Foam::plic::face_from_lbls
(
    const label a,
    const label b,
    const label c,
    const label d
)
{
    face tFc(4);
    tFc[0] = a; tFc[1] = b; tFc[2] = c; tFc[3] = d;
    return tFc;
}


void Foam::plic::makePatchFaceCellInfo
(
    const fvPatch& curPatch,
    const label faceI,
    cellInfo& curCellInfo
)
{
    const polyPatch& pp = curPatch.patch();
    const pointField& patchPoints = pp.localPoints();    
    const face& curFace = pp.localFaces()[faceI];
    vectorField faceCellCtrs(pp.faceCellCentres());

    vector transln_vec = 2.0*(pp.faceCentres()[faceI] - faceCellCtrs[faceI]); 

    pointField& cellPoints = curCellInfo.points();
    faceList& cellFaces = curCellInfo.faces();

    label nCellPts = 0;
    cellPoints.setSize(8);
    forAll(curFace, pointI)
    {
        cellPoints[nCellPts++] = patchPoints[curFace[pointI]];
    }
    forAll(curFace, pointI)
    {
        cellPoints[nCellPts++] = patchPoints[curFace[pointI]] + transln_vec;
    }

    cellFaces.setSize(6);

    face face0(4);
    face0[0] = 0; face0[1] = 1; face0[2] = 2; face0[3] = 3;
    cellFaces[0] = face0;

    face face1(4);
    face1[0] = 1; face1[1] = 0; face1[2] = 4; face1[3] = 5;
    cellFaces[1] = face1;

    face face2(4);
    face2[0] = 2; face2[1] = 1; face2[2] = 5; face2[3] = 6;
    cellFaces[2] = face2;

    face face3(4);
    face3[0] = 3; face3[1] = 2; face3[2] = 6; face3[3] = 7;
    cellFaces[3] = face3;

    face face4(4);
    face4[0] = 0; face4[1] = 3; face4[2] = 7; face4[3] = 4;
    cellFaces[4] = face4;

    face face5(4);
    face5[0] = 5; face5[1] = 4; face5[2] = 7; face5[3] = 6;
    cellFaces[5] = face5;

    curCellInfo.calc_centreAndVol();
    curCellInfo.correct_faceDir();
}


void Foam::plic::cellInfoAndPlanes_collectData()
{
    const pointField& points = mesh().points();
    const faceList& faces = mesh().faces();
    const labelList& own = mesh().faceOwner();
    const cellList& cells = mesh().cells();
    const volVectorField& C = mesh().C();
    const mapDistribute& map = faceStencil().map();

    if(debug2_)
    {
        Info<< "========================================================================" << nl
            << "                      Making cellInfo for cells" << nl
            << "========================================================================" << nl 
            << endl;
    }

    // Make and insert my internal values
    forAll(cells, cellI)
    {
        cellInfo curCellInfo(points, faces, own, cells[cellI], cellI);
        cells_flatFld_[cellI] = curCellInfo;
        meshC_.internalField()[cellI] = C.internalField()[cellI];
        if(debug2_)
        {
            Info<< "cellInfo cell " << cellI << endl;
            Foam::plicFuncs::display_cellInfo(curCellInfo);
        }
    }

    if(debug2_)
    {
        Info<< "========================================================================" << nl
            << "                  Done making cellInfo for cells" << nl
            << "========================================================================" << nl 
            << endl;
    }

    if(debug2_)
    {
        Info<< "========================================================================" << nl
            << "              Making cellInfo for boundary face cells" << nl
            << "========================================================================" << nl 
            << endl;
    }

    // Make and insert my boundary values
    forAll(mesh().boundary(), patchI)
    {
        const fvPatch& pp = mesh().boundary()[patchI];
        fvPatchVectorField& pmeshC = meshC_.boundaryField()[patchI];

        if(debug2_)
        {
            Info<< "--------------------------------------------------------------------" << nl
                << "                          Patch " << patchI << nl
                << "--------------------------------------------------------------------" << nl 
                << endl;
        }
        
        label nCompact =
            pp.start()
            -mesh().nInternalFaces()
            +mesh().nCells();
        
        cellInfo curCellInfo;

        forAll(pp, faceI)
        {
            makePatchFaceCellInfo(pp, faceI, curCellInfo);

            if(debug2_)
            {
                Info<< "patch face index: " << faceI << "    flat fld index: " << nCompact << endl;
                Foam::plicFuncs::display_cellInfo(curCellInfo);
            }

            pmeshC[faceI] = curCellInfo.centre();
            cells_flatFld_[nCompact++] = curCellInfo;
        }

        if(debug2_)
        {
            Info<< "--------------------------------------------------------------------" << nl
                << "                       Done Patch " << patchI << nl
                << "--------------------------------------------------------------------" << nl 
                << endl;
        }        
    }

    if(debug2_)
    {
        Info<< "========================================================================" << nl
            << "           Done making cellInfo for boundary face cells" << nl
            << "========================================================================" << nl 
            << endl;
    }
    
    if(debug2_)
    {
        Info<< "========================================================================" << nl
            << "                Making cell bounding planes list" << nl
            << "========================================================================" << nl 
            << endl;
    }

    Plane initPlane(vector::one);

    //label n_flatFld_local = mesh().nCells() + mesh().nFaces() - mesh().nInternalFaces();

    forAll(cells, i)
    {
        const cellInfo& curCellInfo = cells_flatFld_[i];
        const pointField& cellPts = curCellInfo.points();
        const faceList& cellFcs = curCellInfo.faces();
        List<Plane> cellPlns(curCellInfo.faces().size(), initPlane);
        
        if(debug2_)
        {
            Info<< "--------------------------------------------------------------------" << nl
                << "            Planes for flat fld cell index: " << i << nl
                << "--------------------------------------------------------------------" << nl 
                << endl;
        }

        for(label plnI=0; plnI<cellFcs.size(); plnI++)
        {
            const face& curFc = cellFcs[plnI];            
            const point& pt1 = cellPts[curFc[0]];
            const point& pt2 = cellPts[curFc[1]];
            const point& pt3 = cellPts[curFc[2]];            
            Plane curPln(pt1, pt2, pt3);

            if(debug2_)
            {
                Info<< "Bounding plane no. " << plnI << endl;
                Foam::plicFuncs::display_Plane(curPln);
            }

            cellPlns[plnI] = curPln;
        }

        cellPlns_flatFld_[i] = cellPlns;

        if(debug2_)
        {
            Info<< "--------------------------------------------------------------------" << nl
                << "           Done planes for flat fld cell index: " << i << nl
                << "--------------------------------------------------------------------" << nl 
                << endl;
        }
    }

    // Make and insert my boundary values
    forAll(mesh().boundary(), patchI)
    {
        const fvPatch& pp = mesh().boundary()[patchI];

        if(debug2_)
        {
            Info<< "--------------------------------------------------------------------" << nl
                << "                          Patch " << patchI << nl
                << "--------------------------------------------------------------------" << nl 
                << endl;
        }
        
        label nCompact =
            pp.start()
            -mesh().nInternalFaces()
            +mesh().nCells();

        forAll(pp, i)
        {
            const cellInfo& curCellInfo = cells_flatFld_[nCompact];
            const pointField& cellPts = curCellInfo.points();
            const faceList& cellFcs = curCellInfo.faces();
            List<Plane> cellPlns(curCellInfo.faces().size(), initPlane);
        
            if(debug2_)
            {
                Info<< "--------------------------------------------------------------------" << nl
                    << "            Planes for flat fld cell index: " << nCompact << nl
                    << "--------------------------------------------------------------------" << nl 
                    << endl;
            }

            for(label plnI=0; plnI<cellFcs.size(); plnI++)
            {
                const face& curFc = cellFcs[plnI];            
                const point& pt1 = cellPts[curFc[0]];
                const point& pt2 = cellPts[curFc[1]];
                const point& pt3 = cellPts[curFc[2]];            
                Plane curPln(pt1, pt2, pt3);

                if(debug2_)
                {
                    Info<< "Bounding plane no. " << plnI << endl;
                    Foam::plicFuncs::display_Plane(curPln);
                }

                cellPlns[plnI] = curPln;
            }

            cellPlns_flatFld_[nCompact] = cellPlns;

            if(debug2_)
            {
                Info<< "--------------------------------------------------------------------" << nl
                    << "           Done planes for flat fld cell index: " << nCompact << nl
                    << "--------------------------------------------------------------------" << nl 
                    << endl;
            }
            nCompact++;
        }
    }

    if(debug2_)
    {
        Info<< "========================================================================" << nl
            << "                Done making cell bounding planes list" << nl
            << "========================================================================" << nl 
            << endl;
    }

    // Do all swapping
    //map.distribute(cells_flatFld_);
    map.distribute(cellPlns_flatFld_);
}


void Foam::plic::intfcInfo_collectData()
{
    const mapDistribute& map = faceStencil().map();

    // Insert my internal values
    forAll(alpha_ph1_, cellI)
    {
        alpha_ph1_flatFld_[cellI] = alpha_ph1_[cellI];
        nHat_flatFld_[cellI] = nHat_[cellI];
        gradAlpha1_flatFld_[cellI] = gradAlpha1_[cellI];
        C_intfc_flatFld_[cellI] = C_intfc_[cellI];
        rho1_flatFld_[cellI] = rho1_[cellI];
        rho0_flatFld_[cellI] = rho0_[cellI];
    }
    // Insert my boundary values
    forAll(alpha_ph1_.boundaryField(), patchI)
    {
        const fvPatchScalarField& palpha = alpha_ph1_.boundaryField()[patchI];
        const fvPatchVectorField& pnHat = nHat_.boundaryField()[patchI];
        const fvPatchVectorField& pgradAlpha1 = gradAlpha1_.boundaryField()[patchI];
        const fvPatchVectorField& pCIntfc = C_intfc_.boundaryField()[patchI];
        const fvPatchScalarField& prho1 = rho1_.boundaryField()[patchI];
        const fvPatchScalarField& prho0 = rho0_.boundaryField()[patchI];

        label nCompact =
            palpha.patch().start()
            -alpha_ph1_.mesh().nInternalFaces()
            +alpha_ph1_.mesh().nCells();

        forAll(palpha, i)
        {
            alpha_ph1_flatFld_[nCompact] = palpha[i];
            nHat_flatFld_[nCompact] = pnHat[i];
            gradAlpha1_flatFld_[nCompact] = pgradAlpha1[i];
            C_intfc_flatFld_[nCompact] = pCIntfc[i];
            rho1_flatFld_[nCompact] = prho1[i];
            rho0_flatFld_[nCompact] = prho0[i];
            nCompact++;
        }
    }

    // Do all swapping
    map.distribute(alpha_ph1_flatFld_);
    map.distribute(nHat_flatFld_);
    map.distribute(gradAlpha1_flatFld_);
    map.distribute(C_intfc_flatFld_);
    map.distribute(rho1_flatFld_);
    map.distribute(rho0_flatFld_);

    if(debug2_)
    {
        Foam::plicFuncs::write_flatFld(alpha_ph1_flatFld_, alpha_ph1_);
        Foam::plicFuncs::write_flatFld(nHat_flatFld_, nHat_);
        Foam::plicFuncs::write_flatFld(gradAlpha1_flatFld_, gradAlpha1_);
        Foam::plicFuncs::write_flatFld(C_intfc_flatFld_, C_intfc_);
        Foam::plicFuncs::write_flatFld(rho1_flatFld_, rho1_);
        Foam::plicFuncs::write_flatFld(rho0_flatFld_, rho0_);
    }
}


void Foam::plic::Y_collectData
(
    const PtrList<volScalarField>& Y1,
    const PtrList<volScalarField>& Y0    
)
{
    const mapDistribute& map = faceStencil().map();

    for(label i=0; i<(nSpecies_ - 1); i++)
    {
        const volScalarField& Y1i = Y1[i];
        const volScalarField& Y0i = Y0[i];
        List<scalar>& Y1i_flatFld = Y1_flatFld_[i];
        List<scalar>& Y0i_flatFld = Y0_flatFld_[i];

        // Insert my internal values
        forAll(Y1i, cellI)
        {        
            Y1i_flatFld[cellI] = Y1i[cellI];
            Y0i_flatFld[cellI] = Y0i[cellI];
        }
        // Insert my boundary values
        forAll(Y1i.boundaryField(), patchI)
        {        
            const fvPatchScalarField& pY1i = Y1i.boundaryField()[patchI];        
            const fvPatchScalarField& pY0i = Y0i.boundaryField()[patchI];      
            label nCompact =
                pY1i.patch().start()
                -Y1i.mesh().nInternalFaces()
                +Y1i.mesh().nCells();

            forAll(pY1i, faceI)
            {           
                Y1i_flatFld[nCompact] = pY1i[faceI];
                Y0i_flatFld[nCompact] = pY0i[faceI];
                nCompact++;
            }
        }

        // Do all swapping    
        map.distribute(Y1i_flatFld);
        map.distribute(Y0i_flatFld);
    
        if(debug2_)
        {        
            Foam::plicFuncs::write_flatFld(Y1i_flatFld, Y1i);
            Foam::plicFuncs::write_flatFld(Y0i_flatFld, Y0i);
        }
    }
}


void Foam::plic::c_h_collectData
(
    const PtrList<volScalarField>& c1,
    const PtrList<volScalarField>& c0,
    const volScalarField& h1,
    const volScalarField& h0
)
{
    const mapDistribute& map = faceStencil().map();

    for(label i=0; i<(nSpecies_ - 1); i++)
    {
        const volScalarField& c1i = c1[i];
        const volScalarField& c0i = c0[i];
        List<scalar>& c1i_flatFld = c1_flatFld_[i];
        List<scalar>& c0i_flatFld = c0_flatFld_[i];

        // Insert my internal values
        forAll(c1i, cellI)
        {        
            c1i_flatFld[cellI] = c1i[cellI];
            c0i_flatFld[cellI] = c0i[cellI];
        }
        // Insert my boundary values
        forAll(c1i.boundaryField(), patchI)
        {        
            const fvPatchScalarField& pc1i = c1i.boundaryField()[patchI];        
            const fvPatchScalarField& pc0i = c0i.boundaryField()[patchI];      
            label nCompact =
                pc1i.patch().start()
                -c1i.mesh().nInternalFaces()
                +c1i.mesh().nCells();

            forAll(pc1i, faceI)
            {           
                c1i_flatFld[nCompact] = pc1i[faceI];
                c0i_flatFld[nCompact] = pc0i[faceI];
                nCompact++;
            }
        }

        // Do all swapping    
        map.distribute(c1i_flatFld);
        map.distribute(c0i_flatFld);
    
        if(debug2_)
        {        
            Foam::plicFuncs::write_flatFld(c1i_flatFld, c1i);
            Foam::plicFuncs::write_flatFld(c0i_flatFld, c0i);
        }
    }    

    // Insert my internal values
    forAll(h1, cellI)
    {        
        h1_flatFld_[cellI] = h1[cellI];
        h0_flatFld_[cellI] = h0[cellI];
    }
    // Insert my boundary values
    forAll(h1.boundaryField(), patchI)
    {        
        const fvPatchScalarField& ph1 = h1.boundaryField()[patchI];        
        const fvPatchScalarField& ph0 = h0.boundaryField()[patchI];      
        label nCompact =
            ph1.patch().start()
            -h1.mesh().nInternalFaces()
            +h1.mesh().nCells();

        forAll(ph1, faceI)
        {           
            h1_flatFld_[nCompact] = ph1[faceI];
            h0_flatFld_[nCompact] = ph0[faceI];
            nCompact++;
        }
    }

    // Do all swapping    
    map.distribute(h1_flatFld_);
    map.distribute(h0_flatFld_);
    
    if(debug2_)
    {        
        Foam::plicFuncs::write_flatFld(h1_flatFld_, h1);
        Foam::plicFuncs::write_flatFld(h0_flatFld_, h0);
    }
}


void Foam::plic::c_collectData
(
    const PtrList<volScalarField>& c1,
    const PtrList<volScalarField>& c0    
)
{
    const mapDistribute& map = faceStencil().map();

    for(label i=0; i<(nSpecies_ - 1); i++)
    {
        const volScalarField& c1i = c1[i];
        const volScalarField& c0i = c0[i];
        List<scalar>& c1i_flatFld = c1_flatFld_[i];
        List<scalar>& c0i_flatFld = c0_flatFld_[i];

        // Insert my internal values
        forAll(c1i, cellI)
        {        
            c1i_flatFld[cellI] = c1i[cellI];
            c0i_flatFld[cellI] = c0i[cellI];
        }
        // Insert my boundary values
        forAll(c1i.boundaryField(), patchI)
        {        
            const fvPatchScalarField& pc1i = c1i.boundaryField()[patchI];        
            const fvPatchScalarField& pc0i = c0i.boundaryField()[patchI];      
            label nCompact =
                pc1i.patch().start()
                -c1i.mesh().nInternalFaces()
                +c1i.mesh().nCells();

            forAll(pc1i, faceI)
            {           
                c1i_flatFld[nCompact] = pc1i[faceI];
                c0i_flatFld[nCompact] = pc0i[faceI];
                nCompact++;
            }
        }

        // Do all swapping    
        map.distribute(c1i_flatFld);
        map.distribute(c0i_flatFld);
    
        if(debug2_)
        {        
            Foam::plicFuncs::write_flatFld(c1i_flatFld, c1i);
            Foam::plicFuncs::write_flatFld(c0i_flatFld, c0i);
        }
    }    
}


Foam::scalar Foam::plic::tet_cell_intersect
(
    tetPoints& curTet,
    const label& cellI
)
{        
    scalar alpha_cellI = alpha_ph1_flatFld_[cellI];
    //scalar nHat_cellI_mag = mag(nHat_flatFld_[cellI]);
    scalar gradAlpha1_cellI_mag = mag(gradAlpha1_flatFld_[cellI]);

    if(debugF_)
    {
        Info<< "Cell alpha1:  " << alpha_cellI << nl << endl; 
    }

    if((alpha_cellI > ALPHA_2PH_MIN) && (alpha_cellI < ALPHA_2PH_MAX) && (gradAlpha1_cellI_mag > GRADALPHA_MIN))
    {
        scalar V_alpha1 = 0;
        List<Plane> cellPlns = cellPlns_flatFld_[cellI];
                
        Plane intfcPln(C_intfc_flatFld_[cellI], nHat_flatFld_[cellI]);
        cellPlns.resize(cellPlns.size()+1);
        cellPlns[cellPlns.size()-1] = intfcPln;        

        tet_slice_with_plns(curTet, cellPlns, V_alpha1, 0);

        if(debugF_)
        {
            Info<< "Flux contribution from cell " << cellI << ":  " << V_alpha1 << nl << endl;
        }

        return V_alpha1;
    } 
    else
    {
        scalar V_alpha1 = 0;
        List<Plane> cellPlns = cellPlns_flatFld_[cellI];

        tet_slice_with_plns(curTet, cellPlns, V_alpha1, 0);

        if(debugF_)
        {
            Info<< "Flux contribution from cell " << cellI << ":  " << V_alpha1 << nl << endl;
        }

        return V_alpha1*alpha_cellI;
    }

    /*    if(alpha_cellI < ALPHA_2PH_MIN)
        {
        if(debugF2_)
        {
        Info<< "Flux contribution from cell " << cellI << ":  " << 0 << nl << endl;
        }

        return 0;
        }
        else
        {
        scalar V_alpha1 = 0;
        List<Plane> cellPlns = cellPlns_flatFld_[cellI];

        if((alpha_cellI < ALPHA_2PH_MAX) && (nHat_cellI_mag > GRADALPHA_MIN))
        {
        Plane intfcPln(C_intfc_flatFld_[cellI], nHat_flatFld_[cellI]);
        cellPlns.resize(cellPlns.size()+1);
        cellPlns[cellPlns.size()-1] = intfcPln;
        }

        tet_slice_with_plns(curTet, cellPlns, V_alpha1, 0);

        if(debugF2_)
        {
        Info<< "Flux contribution from cell " << cellI << ":  " << V_alpha1 << nl << endl;
        }

        return V_alpha1;
        }
        */
}


void Foam::plic::tet_cell_intersect
(
    tetPoints& curTet,
    const label& cellI,
    scalar& V_alpha1,
    scalar& V_alpha0
)
{        
    scalar alpha_cellI = alpha_ph1_flatFld_[cellI];
    //scalar nHat_cellI_mag = mag(nHat_flatFld_[cellI]);
    scalar gradAlpha1_cellI_mag = mag(gradAlpha1_flatFld_[cellI]);

    if(debugF2_)
    {
        Info<< "Cell alpha1:  " << alpha_cellI << nl << endl; 
    }

    scalar V_tot = 0;
    V_alpha1 = 0;
    V_alpha0 = 0;

    if((alpha_cellI > ALPHA_2PH_MIN) && (alpha_cellI < ALPHA_2PH_MAX) && (gradAlpha1_cellI_mag > GRADALPHA_MIN))
    {
        List<Plane> cellPlns = cellPlns_flatFld_[cellI];
                
        Plane intfcPln(C_intfc_flatFld_[cellI], nHat_flatFld_[cellI]);
        cellPlns.resize(cellPlns.size()+1);
        cellPlns[cellPlns.size()-1] = intfcPln;        

        tet_slice_with_plns(curTet, cellPlns, V_alpha1, V_tot, 0);

        V_alpha0 = V_tot - V_alpha1;
    } 
    else
    {
        List<Plane> cellPlns = cellPlns_flatFld_[cellI];

        tet_slice_with_plns(curTet, cellPlns, V_tot, 0);        

        V_alpha1 = V_tot*alpha_cellI;
        V_alpha0 = V_tot - V_alpha1;
    }
}


void Foam::plic::tet_slice_with_plns
(
    tetPoints& curTet,
    const List<Plane>& cellPlns,
    scalar& V_alpha1,
    const label& pln_lbl
)
{
    static tetPointRef::tetIntersectionList insideTets;
    label nInside = 0;
    tetPointRef::storeOp inside(insideTets, nInside);
    tetPointRef::dummyOp outside;

    plane curPln(cellPlns[pln_lbl].refPoint(), cellPlns[pln_lbl].normal());

    if(debugF2_)
    {
        Info<< "Slicing tet with plane level: " << pln_lbl << nl
            << "Plane normal: " << curPln.normal() << "  Plane refPoint: " << curPln.refPoint() << nl
            << "Tet to be sliced:" << nl
            << "Points:  " << curTet.tet().a() << "  " << curTet.tet().b() << "  " << curTet.tet().c() << "  " << curTet.tet().d() << nl
            << "Volume:  " << curTet.tet().mag() << endl;
    }

    curTet.tet().sliceWithPlane(curPln, inside, outside);

    if(debugF2_)
    {
        Info<< nl
            << "No. of tets in plane half-space: " << nInside << nl << endl;

        for(label tetI=0; tetI<=nInside; tetI++)
        {
            if(debugF2_)
            {                
                Info<< "Tetrahedron " << tetI << nl
                    << "Points: " << insideTets[tetI].tet().a() << "  " << insideTets[tetI].tet().b() << "  " << insideTets[tetI].tet().c() << "  " << insideTets[tetI].tet().d() << nl
                    << "Volume: " << insideTets[tetI].tet().mag() << nl << endl;
            }
        }
    }

    List<tetPoints> newTets(nInside);
    label nNewTets = 0;

    for(label tetI=0; tetI<nInside; tetI++)
    {
        if(mag(insideTets[tetI].tet().mag()) > SMALLEST_TET_VOL)
        {
            newTets[nNewTets++] = insideTets[tetI];
        }
    }

    newTets.resize(nNewTets);

    for(label tetI=0; tetI<nNewTets; tetI++)
    {
        if(debugF2_)
        {                
            Info<< "New tetrahedron " << tetI << nl
                << "Points: " << newTets[tetI].tet().a() << "  " << newTets[tetI].tet().b() << "  " << newTets[tetI].tet().c() << "  " << newTets[tetI].tet().d() << nl
                << "Volume: " << newTets[tetI].tet().mag() << nl << endl;
        }
    }

    if(nNewTets > 0)
    {
        if(pln_lbl == (cellPlns.size()-1))
        {
            if(debugF2_)
            {
                Info<< "Done slicing tet with last plane (cell bounding plane or intfc plane)..." << nl
                    << "Adding resulting tet volumes to flux contribution from original tet..." << nl << endl;
            }

            for(label tetI=0; tetI<nNewTets; tetI++)
            {                
                V_alpha1 += mag(newTets[tetI].tet().mag());

                if(debugF2_)
                {
                    Info<< "Adding volume of final level tet " << tetI << nl
                        << "Volume added:  " << newTets[tetI].tet().mag() << "    " << mag(newTets[tetI].tet().mag()) << nl
                        << "V_alpha1 (total so far):  " << V_alpha1 << nl << endl;
                }
            }
        }
        else
        {
            for(label tetI=0; tetI<nNewTets; tetI++)
            {
                if(debugF2_)
                {
                    Info<< "Slicing level " << pln_lbl << " produced tet " << tetI << " with level " << pln_lbl+1 << " plane" << nl
                        << "Points: " << newTets[tetI].tet().a() << "  " << newTets[tetI].tet().b() << "  " << newTets[tetI].tet().c() << "  " << newTets[tetI].tet().d() << nl
                        << "Volume: " << newTets[tetI].tet().mag() << nl << endl;
                }

                tet_slice_with_plns(newTets[tetI], cellPlns, V_alpha1, pln_lbl+1);
            }
        }
    }
    else
    {
        if(debugF2_)
        {
            Info<< "Null intersection with sub-cell for this tet..." << nl << endl;
        }
    }
}


void Foam::plic::tet_slice_with_plns
(
    tetPoints& curTet,
    const List<Plane>& cellPlns,
    scalar& V_alpha1,
    scalar& V_tot,
    const label& pln_lbl
)
{
    static tetPointRef::tetIntersectionList insideTets;
    label nInside = 0;
    tetPointRef::storeOp inside(insideTets, nInside);
    tetPointRef::dummyOp outside;

    plane curPln(cellPlns[pln_lbl].refPoint(), cellPlns[pln_lbl].normal());

    if(debugF2_)
    {
        Info<< "Slicing tet with plane level: " << pln_lbl << nl
            << "Plane normal: " << curPln.normal() << "  Plane refPoint: " << curPln.refPoint() << nl
            << "Tet to be sliced:" << nl
            << "Points:  " << curTet.tet().a() << "  " << curTet.tet().b() << "  " << curTet.tet().c() << "  " << curTet.tet().d() << nl
            << "Volume:  " << curTet.tet().mag() << endl;
    }

    curTet.tet().sliceWithPlane(curPln, inside, outside);

    if(debugF2_)
    {
        Info<< nl
            << "No. of tets in plane half-space: " << nInside << nl << endl;

        for(label tetI=0; tetI<=nInside; tetI++)
        {
            if(debugF2_)
            {                
                Info<< "Tetrahedron " << tetI << nl
                    << "Points: " << insideTets[tetI].tet().a() << "  " << insideTets[tetI].tet().b() << "  " << insideTets[tetI].tet().c() << "  " << insideTets[tetI].tet().d() << nl
                    << "Volume: " << insideTets[tetI].tet().mag() << nl << endl;
            }
        }
    }

    List<tetPoints> newTets(nInside);
    label nNewTets = 0;

    for(label tetI=0; tetI<nInside; tetI++)
    {
        if(mag(insideTets[tetI].tet().mag()) > SMALLEST_TET_VOL)
        {
            newTets[nNewTets++] = insideTets[tetI];
        }
    }

    newTets.resize(nNewTets);

    for(label tetI=0; tetI<nNewTets; tetI++)
    {
        if(debugF2_)
        {                
            Info<< "New tetrahedron " << tetI << nl
                << "Points: " << newTets[tetI].tet().a() << "  " << newTets[tetI].tet().b() << "  " << newTets[tetI].tet().c() << "  " << newTets[tetI].tet().d() << nl
                << "Volume: " << newTets[tetI].tet().mag() << nl << endl;
        }
    }

    if(nNewTets > 0)
    {
        if(pln_lbl == (cellPlns.size()-1))
        {
            if(debugF2_)
            {
                Info<< "Done slicing tet with intfc plane..." << nl
                    << "Adding resulting tet volumes to ph-1 flux contribution from given cell..." << nl << endl;
            }

            for(label tetI=0; tetI<nNewTets; tetI++)
            {                
                V_alpha1 += mag(newTets[tetI].tet().mag());

                if(debugF2_)
                {
                    Info<< "Adding volume of final level tet " << tetI << nl
                        << "Volume added:  " << newTets[tetI].tet().mag() << "    " << mag(newTets[tetI].tet().mag()) << nl
                        << "V_alpha1 (total so far):  " << V_alpha1 << nl << endl;
                }
            }
        }
        else if(pln_lbl == (cellPlns.size()-2))
        {
            if(debugF2_)
            {
                Info<< "Done slicing tet with last cell bounding plane..." << nl
                    << "Adding resulting tet volumes to flux contribution from given cell..." << nl << endl;
            }

            for(label tetI=0; tetI<nNewTets; tetI++)
            {                
                V_tot += mag(newTets[tetI].tet().mag());

                if(debugF2_)
                {
                    Info<< "Adding volume of penultimate level tet " << tetI << nl
                        << "Volume added:  " << newTets[tetI].tet().mag() << "    " << mag(newTets[tetI].tet().mag()) << nl
                        << "V_tot (total so far):  " << V_tot << nl << endl;
                }

                tet_slice_with_plns(newTets[tetI], cellPlns, V_alpha1, V_tot, pln_lbl+1);
            }
        }
        else
        {
            for(label tetI=0; tetI<nNewTets; tetI++)
            {
                if(debugF2_)
                {
                    Info<< "Slicing level " << pln_lbl << " produced tet " << tetI << " with level " << pln_lbl+1 << " plane" << nl
                        << "Points: " << newTets[tetI].tet().a() << "  " << newTets[tetI].tet().b() << "  " << newTets[tetI].tet().c() << "  " << newTets[tetI].tet().d() << nl
                        << "Volume: " << newTets[tetI].tet().mag() << nl << endl;
                }

                tet_slice_with_plns(newTets[tetI], cellPlns, V_alpha1, V_tot, pln_lbl+1);
            }
        }
    }
    else
    {
        if(debugF2_)
        {
            Info<< "Null intersection with sub-cell for this tet..." << nl << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::plic::plic
(
    const fvMesh& mesh,
    const volScalarField& alpha_ph1,
    const volScalarField& alpha_ph0,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volScalarField& rho1,
    const volScalarField& rho0,
    const volScalarField& rho,
    OFstream& osPlic
)
    :
    mesh_(mesh),
    faceStencil_(mesh),
    lsGrad_(mesh),
    alpha_ph1_(alpha_ph1),
    alpha_ph0_(alpha_ph0),
    rAlpha0_
    (
        IOobject
        (
            "rAlpha0",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rAlpha0_0", dimless, 1.0)
    ),
    transPropDict_
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    U_(U),
    ptInterp_(mesh),
    U_pts_
    (
        IOobject
        (
            "U_pts",
            mesh.polyMesh::instance(),
            mesh
        ),
        pointMesh::New(mesh),
        dimensionedVector("U_pts", dimVelocity, vector::zero)
    ),
    alpha_ph1_pts_
    (
        IOobject
        (
            "alpha_ph1_pts",
            mesh.polyMesh::instance(),
            mesh
        ),
        pointMesh::New(mesh),
        dimensionedScalar("alpha_ph1_pts", dimless, 0)
    ),
    alpha1f_
    (
        IOobject
        (
            "alpha1f",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(alpha_ph1_)
    ),
    rho1_(rho1),
    rho0_(rho0),
    rho_(rho),
    phi_(phi),
    phiAlpha1_
    (
        IOobject
        (
            "phiAlpha1",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phi
    ),
    phiAlpha0_
    (
        IOobject
        (
            "phiAlpha0",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi
    ),
    nMixCells_(0),
    nMixFcs_(0),
    gradAlpha1_
    (
        IOobject
        (
            "gradAlpha1",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha_ph1.mesh(),
        dimensionedVector("gradAlpha1", dimless/dimLength, vector::one),
        zeroGradientFvPatchField<vector>::typeName
    ),
    nHat_
    (
        IOobject
        (
            "nHat",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha_ph1.mesh(),
        dimensionedVector("nHat", dimless, vector::one),
        zeroGradientFvPatchField<vector>::typeName
    ),  
    nHatf_
    (
        IOobject
        (
            "nHatf_plic",
            alpha_ph1.time().timeName(),
            alpha_ph1.mesh()
        ),
        alpha_ph1.mesh(),
        dimensionedScalar("nHatf_plic", dimArea, 0.0)
    ),
    C_intfc_
    (
        IOobject
        (
            "C_intfc",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    A_intfc_
    (
        IOobject
        (
            "A_intfc",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimArea, 0)
    ),
    C_ph1_
    (
        IOobject
        (
            "C_ph1",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    C_ph0_
    (
        IOobject
        (
            "C_ph0",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    meshC_
    (
        IOobject
        (
            "meshC",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    species_(transPropDict_.lookup("species")),
    nSpecies_(species_.size()),
    fAlpha1_(fvc::interpolate(alpha_ph1_,"alpha1")),
    frho1_(fvc::interpolate(rho1_,"rho")),
    frho0_(fvc::interpolate(rho0_,"rho")),
    fY1_(nSpecies_ - 1),
    fY0_(nSpecies_ - 1),
    fc1_(nSpecies_ - 1),
    fc0_(nSpecies_ - 1),
    fh1_
    (
        IOobject
        (
            "fh1",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zeroh", dimMass/dimLength/dimTime/dimTime, 0)
    ),
    fh0_
    (
        IOobject
        (
            "fh0",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zeroh", dimMass/dimLength/dimTime/dimTime, 0)
    ),
    os(osPlic)
{
    //Foam::plicFuncs::write_stencil(faceStencil().stencil(), mesh, "centredFPCCellToFaceStencil");

    const label nCells = mesh.nCells();
    const label nFaces = mesh.nFaces();

    cell_phaseState_.setSize(nCells);
    face_phaseState_own_.setSize(nFaces);
    face_phaseState_nei_.setSize(nFaces);
    cell_near_intfc_.setSize(nCells);
    fc_2ph_flux_needed_.setSize(nFaces);
    Cf_ph1_own_.setSize(nFaces);
    Cf_ph0_own_.setSize(nFaces);
    Af_ph1_own_.setSize(nFaces);
    Af_ph0_own_.setSize(nFaces);
    Cf_ph1_nei_.setSize(nFaces);
    Cf_ph0_nei_.setSize(nFaces);
    Af_ph1_nei_.setSize(nFaces);
    Af_ph0_nei_.setSize(nFaces);

    const surfaceVectorField& meshCf = mesh.Cf();
    const surfaceVectorField& meshSf = mesh.Sf();
    const surfaceScalarField& meshMagSf = mesh.magSf();
    Cf_.setSize(nFaces);
    Sf_.setSize(nFaces);
    magSf_.setSize(nFaces);

    for(label faceI=0; faceI<mesh.nInternalFaces(); faceI++)
    {
        Cf_[faceI] = meshCf[faceI];
        Sf_[faceI] = meshSf[faceI];        
        magSf_[faceI] = meshMagSf[faceI];        
    }
    
    forAll(meshSf.boundaryField(), patchI)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchI];
        const fvsPatchVectorField& pCf = meshCf.boundaryField()[patchI];
        const fvsPatchVectorField& pSf = meshSf.boundaryField()[patchI];
        const fvsPatchScalarField& pMagSf = meshMagSf.boundaryField()[patchI];
        label faceI = pp.start();

        forAll(pSf, fcI)
        {
            Cf_[faceI] = pCf[fcI];
            Sf_[faceI] = pSf[fcI];        
            magSf_[faceI] = pMagSf[fcI];
            faceI++;
        }
    }

    const label nflatFld = faceStencil().map().constructSize();
    cells_flatFld_.setSize(nflatFld);
    cellPlns_flatFld_.setSize(nflatFld);
    alpha_ph1_flatFld_.setSize(nflatFld);
    nHat_flatFld_.setSize(nflatFld);
    gradAlpha1_flatFld_.setSize(nflatFld);
    C_intfc_flatFld_.setSize(nflatFld);

    forAll(fY1_, i)
    {
        fY1_.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "fY1"+Foam::name(i),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zeroY", dimless, 0)
            )
        );

        fY0_.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "fY0"+Foam::name(i),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zeroY", dimless, 0)
            )
        );

        fc1_.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "fc1"+Foam::name(i),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zeroc", dimDensity, 0)
            )
        );

        fc0_.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "fc0"+Foam::name(i),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zeroc", dimDensity, 0)
            )
        );
    }

    rho1_flatFld_.setSize(nflatFld);
    rho0_flatFld_.setSize(nflatFld);
    h1_flatFld_.setSize(nflatFld);
    h0_flatFld_.setSize(nflatFld);
    rAlpha0_flatFld_.setSize(nflatFld);
    
    Y1_flatFld_.setSize(nSpecies_ - 1);
    Y0_flatFld_.setSize(nSpecies_ - 1);
    c1_flatFld_.setSize(nSpecies_ - 1);
    c0_flatFld_.setSize(nSpecies_ - 1);

    for(label i=0; i<(nSpecies_ - 1); i++)
    {
        Y1_flatFld_[i].setSize(nflatFld);
        Y0_flatFld_[i].setSize(nflatFld);
        c1_flatFld_[i].setSize(nflatFld);
        c0_flatFld_[i].setSize(nflatFld);
    }

    ALPHA_2PH_MIN = transPropDict_.lookupOrDefault("ALPHA_2PH_MIN", 1E-09);
    ALPHA_2PH_MAX = 1 - ALPHA_2PH_MIN;
    GRADALPHA_MIN = transPropDict_.lookupOrDefault("GRADALPHA_MIN", 1E-03);
    ALPHA_TOL = transPropDict_.lookupOrDefault("ALPHA_TOL", 1E-12);
    DIST_TOL = transPropDict_.lookupOrDefault("DIST_TOL", 1E-15);
    PT_DIST_TOL = transPropDict_.lookupOrDefault("PT_DIST_TOL", 1E-15);
    SMALLEST_TET_VOL = transPropDict_.lookupOrDefault("SMALLEST_TET_VOL", 1E-24);
    SMALLEST_PHI = transPropDict_.lookupOrDefault("SMALLEST_PHI", 1E-24);
    MAX_BRENT_ITERS = transPropDict_.lookupOrDefault("MAX_BRENT_ITERS", 108);

    Info<< "Tolerance values for plic:" << nl
        << "ALPHA_2PH_MIN:    " << ALPHA_2PH_MIN << nl
        << "GRADALPHA_MIN:    " << GRADALPHA_MIN << nl
        << "ALPHA_TOL:        " << ALPHA_TOL << nl
        << "DIST_TOL:         " << DIST_TOL << nl
        << "PT_DIST_TOL:      " << PT_DIST_TOL << nl 
        << "MAX_BRENT_ITERS:  " << MAX_BRENT_ITERS << nl << endl;

    debug_ = transPropDict_.lookupOrDefault("plic_debug", false);
    debug2_ = transPropDict_.lookupOrDefault("plic_debug2", false);
    debugIR_ = transPropDict_.lookupOrDefault("plic_debugIR", false);
    debugIR2_ = transPropDict_.lookupOrDefault("plic_debugIR2", false);
    debugF_ = transPropDict_.lookupOrDefault("plic_debugF", false);
    debugF2_ = transPropDict_.lookupOrDefault("plic_debugF2", false);    

    brent_iters_max_ = 0;

    if(debug_)
    {
        Info<< "//====================================================================\\" << nl
            << "             Begin cellInfo and Planes construction" << nl
            << "\\====================================================================//" << nl << endl;
    }

    cellInfoAndPlanes_collectData();

    //lsGrad_.correctLeastSquaresVectors(meshC_);

    if(debug_)
    {
        Info<< "//====================================================================\\" << nl
            << "             Done cellInfo and Planes construction" << nl
            << "\\====================================================================//" << nl << endl;
    }

    if(debug_)
    {
        Info<< "//====================================================================\\" << nl
            << "                      Begin Interface correct" << nl
            << "\\====================================================================//" << nl << endl;
    }

    intfc_correct();

    if(debug_)
    {
        Info<< "//====================================================================\\" << nl
            << "                      Done Interface correct" << nl
            << "\\====================================================================//" << nl << endl;
    }

    /*    if(debug_)
        {
        Info<< "//====================================================================\\" << nl
        << "                   Begin 2-phase flux calculation" << nl
        << "\\====================================================================//" << nl << endl;
        }

        calc_face_phaseFluxes();

        if(debug_)
        {
        Info<< "//====================================================================\\" << nl
        << "                   Done 2-phase flux calculation" << nl
        << "\\====================================================================//" << nl << endl;
        }
        */
}


// * * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * //

Foam::plic::~plic()
{}


// * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * //

void Foam::plic::intfc_normal_correct_lsq()
{
    volVectorField gradAlpha(lsGrad_.calcGrad(alpha_ph1_, "gradAlpha1"));

    if(debug2_)
    {
        Foam::plicFuncs::write_field(gradAlpha);
    }    

    //nHat_ = gradAlpha/(mag(gradAlpha) + deltaN_);    

    forAll(nHat_.internalField(), cellI)
    {
        scalar magGradAlpha = mag(gradAlpha.internalField()[cellI]);
        if(magGradAlpha < VSMALL)
        {
            magGradAlpha += VSMALL;
        }
        nHat_.internalField()[cellI] = gradAlpha.internalField()[cellI]/magGradAlpha;
    }

    forAll(nHat_.boundaryField(), patchI)
    {
        const fvPatchVectorField& pgradAlpha = gradAlpha.boundaryField()[patchI];
        fvPatchVectorField& pnHat = nHat_.boundaryField()[patchI];

        forAll(pnHat, faceI)
        {
            scalar magGradAlpha = mag(pgradAlpha[faceI]);
            if(magGradAlpha < VSMALL)
            {
                magGradAlpha += VSMALL;
            }
            pnHat[faceI] = pgradAlpha[faceI]/magGradAlpha;
        }
    }    

    if(debug2_)
    {
        Foam::plicFuncs::write_field(nHat_);
    }
}


void Foam::plic::intfc_normal_correct()
{
    const polyBoundaryMesh& patches = mesh().boundaryMesh();
    
    if(debug_)
    {
        Info<< "Begin alpha1 interpolation to mesh points..." << endl;
    }

    //alpha_ph1_pts_ = ptInterp_.interpolate(alpha_ph1_);
    ptInterp_.interpolate(alpha_ph1_, alpha_ph1_pts_);

    if(debug_)
    {
        Info<< "Done alpha1 interpolation to mesh points..." << endl;
    }

    
    const pointField& pts = mesh_.points();
    const faceList& fcs = mesh_.faces();

    for(label faceI=0; faceI<mesh_.nInternalFaces(); faceI++)
    {
        const labelList& curFc = fcs[faceI];
        label nPoints = curFc.size();

        point fCtr = pts[curFc[0]];
        scalar alpha_fCtr = alpha_ph1_pts_[curFc[0]];
        for (label i=1; i < nPoints; i++)
        {
            fCtr += pts[curFc[i]];
            alpha_fCtr += alpha_ph1_pts_[curFc[i]];
        }
        fCtr /= nPoints;
        alpha_fCtr /= nPoints;

        scalar sumA = 0.0;
        scalar sumAlphaA = 0.0;

        for (label i=0; i < nPoints; i++)
        {
            label curPt_lbl = curFc[i];
            label nextPt_lbl = curFc[(i + 1) % nPoints];
            const point& curPt = pts[curPt_lbl];
            const point& nextPt = pts[nextPt_lbl];

            vector n = (nextPt - curPt)^(fCtr - curPt);
            scalar a = mag(n);
            scalar talpha = (alpha_ph1_pts_[curPt_lbl] + alpha_ph1_pts_[nextPt_lbl] + alpha_fCtr)/3.0;

            sumA += a;
            sumAlphaA += a*talpha;
        }

        alpha1f_[faceI] = sumAlphaA/sumA;
    }

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        fvsPatchScalarField& pAlpha1f = alpha1f_.boundaryField()[patchI];
        label faceI = pp.start();

        forAll(pAlpha1f, fi)
        {
            const labelList& curFc = fcs[faceI];
            label nPoints = curFc.size();

            point fCtr = pts[curFc[0]];
            scalar alpha_fCtr = alpha_ph1_pts_[curFc[0]];
            for (label i=1; i < nPoints; i++)
            {
                fCtr += pts[curFc[i]];
                alpha_fCtr += alpha_ph1_pts_[curFc[i]];
            }
            fCtr /= nPoints;
            alpha_fCtr /= nPoints;

            scalar sumA = 0.0;
            scalar sumAlphaA = 0.0;

            for (label i=0; i < nPoints; i++)
            {
                label curPt_lbl = curFc[i];
                label nextPt_lbl = curFc[(i + 1) % nPoints];
                const point& curPt = pts[curPt_lbl];
                const point& nextPt = pts[nextPt_lbl];

                vector n = (nextPt - curPt)^(fCtr - curPt);
                scalar a = mag(n);
                scalar talpha = (alpha_ph1_pts_[curPt_lbl] + alpha_ph1_pts_[nextPt_lbl] + alpha_fCtr)/3.0;

                sumA += a;
                sumAlphaA += a*talpha;
            }

            pAlpha1f[fi] = sumAlphaA/sumA;

            faceI++;
        }
    }

    if(debug_)
    {
        Info<< "Done alpha1f calculation..." << endl;
    }

    gradAlpha1_ = fvc::grad(alpha1f_);
    gradAlpha1_.correctBoundaryConditions();

    if(debug_)
    {
        Info<< "Done gradAlpha1 calculation..." << endl;
    }

    if(debug2_)
    {
        Foam::plicFuncs::write_field(gradAlpha1_);
    }    
    //Foam::plicFuncs::write_field_vector(gradAlpha1_);

    forAll(nHat_.internalField(), cellI)
    {
        scalar magGradAlpha1 = mag(gradAlpha1_.internalField()[cellI]);
        if(magGradAlpha1 < VSMALL)
        {
            nHat_.internalField()[cellI] = gradAlpha1_.internalField()[cellI]/(magGradAlpha1 + VSMALL);            
        }
        else
        {
            nHat_.internalField()[cellI] = gradAlpha1_.internalField()[cellI]/magGradAlpha1;
        }
    }

    forAll(nHat_.boundaryField(), patchI)
    {
        const fvPatchVectorField& pgradAlpha1 = gradAlpha1_.boundaryField()[patchI];
        fvPatchVectorField& pnHat = nHat_.boundaryField()[patchI];

        forAll(pnHat, faceI)
        {
            scalar magGradAlpha1 = mag(pgradAlpha1[faceI]);
            if(magGradAlpha1 < VSMALL)
            {
                pnHat[faceI] = pgradAlpha1[faceI]/(magGradAlpha1 + VSMALL);
            }
            else
            {
                pnHat[faceI] = pgradAlpha1[faceI]/magGradAlpha1;
            }
        }
    } 

    //nHat_ = gradAlpha1_/(mag(gradAlpha1_) + deltaN_);        
    nHat_.correctBoundaryConditions();    

    nHatf_ = fvc::interpolate(nHat_) & mesh().Sf();

    if(debug2_)
    {
        Foam::plicFuncs::write_field(nHat_);
    }
    //Foam::plicFuncs::write_field_vector(nHat_);
}


void Foam::plic::intfc_correct()
{
    brent_iters_max_ = 0;
    brent_err_max_ = 0;
    brent_max_cell_ = 0;
    brent_max_cellC_ = vector::zero;
    bool WARN_LOW_GRAD = 0;

    if(debug_)
    {
        os<< "Correcting interface normal" << endl;
    }

    intfc_normal_correct();

    if(debug_)
    {
        os<< "Done correcting interface normal" << endl;
    }

    const labelList& own = mesh().owner();
    const labelList& nei = mesh().neighbour();

    nMixCells_ = 0; 

    if(debugIR_)
    {
        os<< "//====================================================================\\" << endl
            << "                           Tagging mixed cells" << endl
            << "\\====================================================================//" << endl << endl
            << "------------------------------------------------------------------------" << endl
            << "cell index        alpha1 value        phase state" << endl
            << "------------------------------------------------------------------------" << endl;
    }
 
    for(label cellI=0; cellI<mesh_.nCells(); cellI++)
    {
        scalar talpha = alpha_ph1_.internalField()[cellI];
        if(talpha < ALPHA_2PH_MIN)
        {
            cell_phaseState_[cellI] = 0;
        }
        else if(talpha > ALPHA_2PH_MAX)
        {
            cell_phaseState_[cellI] = 1;
        }
        else
        {
            if(mag(gradAlpha1_.internalField()[cellI]) > GRADALPHA_MIN) 
            {
                cell_phaseState_[cellI] = 2;
                nMixCells_++;            
            }     
            else
            {
                if(talpha > 0.5)
                {
                    cell_phaseState_[cellI] = 1;
                }
                else
                {
                    cell_phaseState_[cellI] = 0;
                }
            }
        }
        
        if(debugIR_)
        {
            os<< "    " << cellI << "              " << talpha << "              " << cell_phaseState_[cellI] << endl;
        }
    }    

    if(debugIR_)
    {
        os<< endl
            << "//====================================================================\\" << endl
            << "                         Done tagging mixed cells" << endl
            << "\\====================================================================//" << endl << endl << endl;        
    }

    if(debugIR_)
    {
        os<< "//====================================================================\\" << endl
            << "                  Interface plane reconstruction in cells" << endl
            << "\\====================================================================//" << endl << endl;        
    }

    nMixFcs_ = 0;
    boolList fc_set(mesh().nFaces(), false);

    const cellList& cells = mesh().cells();
    const faceList& faces = mesh().faces();
    const pointField& points = mesh().points();
    const vectorField& meshC = mesh().C().internalField();
    const surfaceVectorField& meshCf = mesh().Cf();
    const surfaceScalarField& meshMagSf = mesh().magSf();

    const scalarField& alpha1Cells = alpha_ph1_.internalField();
    const vectorField& nHatCells = nHat_.internalField();
    vectorField& C_intfcCells = C_intfc_.internalField();
    scalarField& A_intfcCells = A_intfc_.internalField();
    vectorField& C_ph1Cells = C_ph1_.internalField();
    vectorField& C_ph0Cells = C_ph0_.internalField();    

    forAll(alpha1Cells, cellI)
    {             
        if(debugIR_)
        {
            os<< "================================================================" << nl
                << "                Reconstructing Intfc in cell " << cellI << nl
                << "================================================================" << nl 
                << endl;
        }

        const cell& curCell = cells[cellI];

        if(cell_phaseState_[cellI] == 2)
        {                        
            cell ph1_cell(0);
            cell ph0_cell(0);
            faceList fcs(0);
            pointField pts(0);
            face intfc_face(0);
            labelList ph1_fcLbls(curCell.size());
            labelList ph0_fcLbls(curCell.size());

            if(mag(nHatCells[cellI]) < GRADALPHA_MIN && WARN_LOW_GRAD == 1)
            {
                os<< "Warning: Interface normal has low magnitude" << nl
                    << "Cell: " << cellI << " alpha1: " << alpha1Cells[cellI] << " nHat: " << nHatCells[cellI] << endl;
            }

            intfc_cell_reconstruct
            (
                nHatCells[cellI],
                alpha1Cells[cellI],
                curCell,
                faces,
                points,
                ph1_cell,
                ph0_cell,
                fcs,
                pts,
                intfc_face,
                ph1_fcLbls,
                ph0_fcLbls
            ); 

            if(brent_iters_tmp_ > brent_iters_max_)
            {
                brent_iters_max_ = brent_iters_tmp_;
                brent_err_max_ = brent_err_tmp_;
                brent_max_cell_ = cellI;
                brent_max_cellC_ = meshC[cellI];
                brent_max_cellAlpha_ = alpha1Cells[cellI];
                brent_max_cell_nHat_ = nHatCells[cellI];
                brent_max_cell_gradAlpha1_ = gradAlpha1_.internalField()[cellI];
            }

            if(intfc_face.size())
            {
                C_intfcCells[cellI] = intfc_face.centre(pts);
                A_intfcCells[cellI] = intfc_face.mag(pts);
            }
            else
            {
                C_intfcCells[cellI] = meshC[cellI];
                A_intfcCells[cellI] = 0;
            }

            C_ph1Cells[cellI] = ph1_cell.centre(pts, fcs);
            C_ph0Cells[cellI] = ph0_cell.centre(pts, fcs);

            if(debugIR_)
            {
                os<< "Interface centre:      " << C_intfcCells[cellI] << endl
                    << "Interface area:        " << A_intfcCells[cellI] << endl
                    << "Phase 1 centroid:      " << C_ph1Cells[cellI] << endl
                    << "Phase 0 centroid:      " << C_ph0Cells[cellI] << endl << endl
                    << "----------------------------------------------------------------" << endl
                    << "                   Phase values on cell faces" << endl
                    << "----------------------------------------------------------------" << endl << endl;
            }

            for(label faceI=0; faceI<curCell.size(); faceI++)
            {
                label curFc_lbl = curCell[faceI];
                label curFcOwn = own[curFc_lbl];

                if(curFcOwn == cellI)
                {
                    if(ph1_fcLbls[faceI] == -1)
                    {
                        face ph0_curFc = fcs[ph0_fcLbls[faceI]];
                        Cf_ph0_own_[curFc_lbl] = ph0_curFc.centre(pts);
                        Af_ph0_own_[curFc_lbl] = ph0_curFc.mag(pts);
                        Cf_ph1_own_[curFc_lbl] = Cf_ph0_own_[curFc_lbl];
                        Af_ph1_own_[curFc_lbl] = 0;
                        face_phaseState_own_[curFc_lbl] = 0;                        
                    }//end if(ph1_fcLbls[faceI] == -1) 
                    else if(ph0_fcLbls[faceI] == -1)
                    {
                        face ph1_curFc = fcs[ph1_fcLbls[faceI]];
                        Cf_ph1_own_[curFc_lbl] = ph1_curFc.centre(pts);
                        Af_ph1_own_[curFc_lbl] = ph1_curFc.mag(pts);
                        Cf_ph0_own_[curFc_lbl] = Cf_ph1_own_[curFc_lbl];
                        Af_ph0_own_[curFc_lbl] = 0;
                        face_phaseState_own_[curFc_lbl] = 1;                        
                    }//end if(ph1_fcLbls[faceI] == -1)
                    else if(ph1_fcLbls[faceI] == -2)
                    {
                        Cf_ph1_own_[curFc_lbl] = faces[curFc_lbl].centre(points);
                        Af_ph1_own_[curFc_lbl] = 0;
                        Cf_ph0_own_[curFc_lbl] = faces[curFc_lbl].centre(points);
                        Af_ph0_own_[curFc_lbl] = 0;
                        face_phaseState_own_[curFc_lbl] = 3;                        
                    }//end if(ph1_fcLbls[faceI] == -1)
                    else
                    {
                        face ph0_curFc = fcs[ph0_fcLbls[faceI]];
                        face ph1_curFc = fcs[ph1_fcLbls[faceI]];                        
                        Cf_ph0_own_[curFc_lbl] = ph0_curFc.centre(pts);
                        Af_ph0_own_[curFc_lbl] = ph0_curFc.mag(pts);
                        Cf_ph1_own_[curFc_lbl] = ph1_curFc.centre(pts);
                        Af_ph1_own_[curFc_lbl] = ph1_curFc.mag(pts);                        
                        face_phaseState_own_[curFc_lbl] = 2;
                    }//end if(ph1_fcLbls[faceI] == -1)
                    if(debugIR_)
                    {
                        os<< "------------------------------------------------------------" << nl
                            << "  Cell face index: " << faceI << "    Global face index: " << curFc_lbl << nl                            
                            << "------------------------------------------------------------" << nl
                            << "Face own: " << curFcOwn << "  cur cell: " << cellI << nl
                            << "Phase 1 centroid: " << Cf_ph1_own_[curFc_lbl] << "    Phase 0 centroid: " << Cf_ph0_own_[curFc_lbl] << nl   
                            << "Phase 1 area: " << Af_ph1_own_[curFc_lbl] << "    Phase 0 area: " << Af_ph0_own_[curFc_lbl] << nl
                            << "Face phase state:  " << face_phaseState_own_[curFc_lbl] << nl
                            << "------------------------------------------------------------" << nl 
                            << endl;
                    }
                }// end if(curFcOwn == cellI)
                else
                {
                    if(ph1_fcLbls[faceI] == -1)
                    {
                        face ph0_curFc = fcs[ph0_fcLbls[faceI]];
                        Cf_ph0_nei_[curFc_lbl] = ph0_curFc.centre(pts);
                        Af_ph0_nei_[curFc_lbl] = ph0_curFc.mag(pts);
                        Cf_ph1_nei_[curFc_lbl] = Cf_ph0_nei_[curFc_lbl];
                        Af_ph1_nei_[curFc_lbl] = 0;
                        face_phaseState_nei_[curFc_lbl] = 0;                        
                    }//end if(ph1_fcLbls[faceI] == -1) 
                    else if(ph0_fcLbls[faceI] == -1)
                    {
                        face ph1_curFc = fcs[ph1_fcLbls[faceI]];
                        Cf_ph1_nei_[curFc_lbl] = ph1_curFc.centre(pts);
                        Af_ph1_nei_[curFc_lbl] = ph1_curFc.mag(pts);
                        Cf_ph0_nei_[curFc_lbl] = Cf_ph1_nei_[curFc_lbl];
                        Af_ph0_nei_[curFc_lbl] = 0;
                        face_phaseState_nei_[curFc_lbl] = 1;                        
                    }//end if(ph1_fcLbls[faceI] == -1)
                    else if(ph1_fcLbls[faceI] == -2)
                    {
                        Cf_ph1_nei_[curFc_lbl] = faces[curFc_lbl].centre(points);
                        Af_ph1_nei_[curFc_lbl] = 0;
                        Cf_ph0_nei_[curFc_lbl] = faces[curFc_lbl].centre(points);
                        Af_ph0_nei_[curFc_lbl] = 0;
                        face_phaseState_nei_[curFc_lbl] = 3;                        
                    }//end if(ph1_fcLbls[faceI] == -1)
                    else
                    {
                        face ph0_curFc = fcs[ph0_fcLbls[faceI]];
                        face ph1_curFc = fcs[ph1_fcLbls[faceI]];                        
                        Cf_ph0_nei_[curFc_lbl] = ph0_curFc.centre(pts);
                        Af_ph0_nei_[curFc_lbl] = ph0_curFc.mag(pts);
                        Cf_ph1_nei_[curFc_lbl] = ph1_curFc.centre(pts);
                        Af_ph1_nei_[curFc_lbl] = ph1_curFc.mag(pts);                        
                        face_phaseState_nei_[curFc_lbl] = 2;
                    }//end if(ph1_fcLbls[faceI] == -1)
                    if(debugIR_)
                    {
                        os<< "------------------------------------------------------------" << nl
                            << "  Cell face index: " << faceI << "    Global face index: " << curFc_lbl << nl                            
                            << "------------------------------------------------------------" << nl
                            << "Face own: " << curFcOwn << "  cur cell: " << cellI << nl
                            << "Phase 1 centroid: " << Cf_ph1_nei_[curFc_lbl] << "    Phase 0 centroid: " << Cf_ph0_nei_[curFc_lbl] << nl   
                            << "Phase 1 area: " << Af_ph1_nei_[curFc_lbl] << "    Phase 0 area: " << Af_ph0_nei_[curFc_lbl] << nl
                            << "Face phase state:  " << face_phaseState_nei_[curFc_lbl] << nl
                            << "------------------------------------------------------------" << nl 
                            << endl;
                    }
                }// end if(curFcOwn == cellI)
            }//end for(label faceI=0; faceI<curCell.size(); faceI++)
        }//end if(cell_phaseState_[cellI] == 2)
        else
        {
            C_intfcCells[cellI] = meshC[cellI];            
            A_intfcCells[cellI] = 0;            
            C_ph1Cells[cellI] = meshC[cellI];
            C_ph0Cells[cellI] = meshC[cellI];

            for(label faceI=0; faceI<curCell.size(); faceI++)
            {
                label curFc_lbl = curCell[faceI];
                label curFcOwn = own[curFc_lbl];

                if(curFcOwn == cellI)
                {
                    Cf_ph0_own_[curFc_lbl] = Cf_[curFc_lbl];
                    Cf_ph1_own_[curFc_lbl] = Cf_[curFc_lbl];
                    if(cell_phaseState_[cellI] == 0)
                    {
                        Af_ph0_own_[curFc_lbl] = magSf_[curFc_lbl];
                        Af_ph1_own_[curFc_lbl] = 0;
                        face_phaseState_own_[curFc_lbl] = 0;
                    }//end if(cell_phaseState_[cellI] == 0)
                    else
                    {
                        Af_ph0_own_[curFc_lbl] = 0;
                        Af_ph1_own_[curFc_lbl] = magSf_[curFc_lbl];
                        face_phaseState_own_[curFc_lbl] = 1;
                    }//end if(cell_phaseState_[cellI] == 0)
                }// end if(curFcOwn == cellI)
                else
                {
                    Cf_ph0_nei_[curFc_lbl] = Cf_[curFc_lbl];
                    Cf_ph1_nei_[curFc_lbl] = Cf_[curFc_lbl];
                    if(cell_phaseState_[cellI] == 0)
                    {
                        Af_ph0_nei_[curFc_lbl] = magSf_[curFc_lbl];
                        Af_ph1_nei_[curFc_lbl] = 0;
                        face_phaseState_nei_[curFc_lbl] = 0;
                    }//end if(cell_phaseState_[cellI] == 0)
                    else
                    {
                        Af_ph0_nei_[curFc_lbl] = 0;
                        Af_ph1_nei_[curFc_lbl] = magSf_[curFc_lbl];
                        face_phaseState_nei_[curFc_lbl] = 1;
                    }//end if(cell_phaseState_[cellI] == 0)
                }// end if(curFcOwn == cellI)
            }//end for(label faceI=0; faceI<curCell.size(); faceI++)
        }//end if(cell_phaseState_[cellI] == 2)        

        if(debugIR_)
        {
            os<< "================================================================" << nl
                << "               Done reconstructing in cell " << cellI << nl
                << "================================================================" << nl
                << endl;
        }
    }//end forAll(alpha1Cells, cellI)    

    if(debugIR_)
    {
        os<< "//====================================================================\\" << nl
            << "                Done interface plane reconstruction in cells" << nl
            << "\\====================================================================//" << nl
            << endl;
    }

    if(debugIR_)
    {
        os<< " Face        " << "Own         " << "alpha1Own        " << "nei          " << "alpha1Nei        " << "phaseStateOwn" << "        phaseStateNei" << endl;
    }

    for(label faceI=0; faceI<mesh().nInternalFaces(); faceI++)
    {
        if((alpha1Cells[own[faceI]] < ALPHA_2PH_MIN && alpha1Cells[nei[faceI]] > ALPHA_2PH_MAX) || (alpha1Cells[nei[faceI]] < ALPHA_2PH_MIN && alpha1Cells[own[faceI]] > ALPHA_2PH_MAX))
        {
            Af_ph0_own_[faceI] = 0;
            Af_ph1_own_[faceI] = 0;
            face_phaseState_own_[faceI] = 3;
            Af_ph0_nei_[faceI] = 0;
            Af_ph1_nei_[faceI] = 0;
            face_phaseState_nei_[faceI] = 3;
        }
        
        if(debugIR_)
        {
            os<< "  " << faceI << "         " << own[faceI] << "      " << alpha1Cells[own[faceI]] << "          " << nei[faceI] << "      " << alpha1Cells[nei[faceI]] << "        " << face_phaseState_own_[faceI] << face_phaseState_nei_[faceI] << endl;
        }
    }

    if(debugIR_)
    {
        os<< "//====================================================================\\" << nl
            << "            Interface plane reconstruction in patch face cells" << nl
            << "\\====================================================================//" << nl 
            << endl;
    }

    const polyBoundaryMesh& patches = mesh().boundaryMesh();
    wordList patchNames(patches.names());

    label nBnd = mesh().nFaces() - mesh().nInternalFaces();

    List<vector> Cf_ph1_bndNei(nBnd);
    List<vector> Cf_ph0_bndNei(nBnd);
    List<scalar> Af_ph1_bndNei(nBnd);
    List<scalar> Af_ph0_bndNei(nBnd);
    List<label> face_phaseState_bndNei(nBnd);

    forAll(alpha_ph1_.boundaryField(), patchI)
    {
        if(debugIR_)
        {
            os<< "=======================================================" << nl
                << "    Patch name:  " << patchNames[patchI] << nl
                << "=======================================================" << nl
                << endl;
        }            

        //May all your coding and simulations work perfectly to bring you Success in all your work.
        const fvPatchScalarField& pAlpha1 = alpha_ph1_.boundaryField()[patchI];
        const fvPatchVectorField& pGradAlpha1 = gradAlpha1_.boundaryField()[patchI];
        const fvPatchVectorField& pnHat = nHat_.boundaryField()[patchI];
        const fvsPatchVectorField& pMeshCf = meshCf.boundaryField()[patchI];
        const fvsPatchScalarField& pMeshMagSf = meshMagSf.boundaryField()[patchI];
        fvPatchVectorField& pC_intfc = C_intfc_.boundaryField()[patchI];
        fvPatchScalarField& pA_intfc = A_intfc_.boundaryField()[patchI];
        fvPatchVectorField& pC_ph1 = C_ph1_.boundaryField()[patchI];
        fvPatchVectorField& pC_ph0 = C_ph0_.boundaryField()[patchI];

        label nfCompact = pAlpha1.patch().start();
        label nCompact =
            nfCompact
            -mesh().nInternalFaces()
            +mesh().nCells();        

        if(isA<zeroGradientFvPatchScalarField>(pAlpha1))
        {
            forAll(pAlpha1, faceI)
            {
                const cellInfo& curCellInfo = cells_flatFld_[nCompact]; 
                cell curCell(6);
                for(label fcI=0; fcI<curCell.size(); fcI++)
                {
                    curCell[fcI] = fcI;
                } 

                if(debugIR_)
                {
                    os<< "=======================================================" << nl
                        << "  Face index: " << nfCompact << "    Flat fld face index: " << nCompact << nl
                        << "=======================================================" << nl
                        << "Boundary alpha1 value:  " << pAlpha1[faceI] << nl                        
                        << "Patch face boundary cell: " << endl;
                    Foam::plicFuncs::display_cell(curCell, curCellInfo.faces(), curCellInfo.points());                        
                }

                if((pAlpha1[faceI] < ALPHA_2PH_MAX) && (pAlpha1[faceI] > ALPHA_2PH_MIN) && (mag(pGradAlpha1[faceI]) > GRADALPHA_MIN))
                {                                                           
                    cell ph1_cell(0);
                    cell ph0_cell(0);
                    faceList fcs(0);
                    pointField pts(0);
                    face intfc_face(0);
                    labelList ph1_fcLbls(curCell.size());
                    labelList ph0_fcLbls(curCell.size());

                    if(mag(pnHat[faceI]) < GRADALPHA_MIN && WARN_LOW_GRAD == 1)
                    {
                        os<< "Warning: Interface normal has low magnitude" << nl
                            << "patch: " << patchI << " face: " << faceI << " alpha1: " << pAlpha1[faceI] << " nHat: " << pnHat[faceI] << endl;
                    }

                    intfc_cell_reconstruct
                    (
                        pnHat[faceI],
                        pAlpha1[faceI],
                        curCell,
                        curCellInfo.faces(),
                        curCellInfo.points(),
                        ph1_cell,
                        ph0_cell,
                        fcs,
                        pts,
                        intfc_face,
                        ph1_fcLbls,
                        ph0_fcLbls
                    );

                    if(brent_iters_tmp_ > brent_iters_max_)
                    {
                        brent_iters_max_ = brent_iters_tmp_;
                        brent_err_max_ = brent_err_tmp_;
                        brent_max_cell_ = nCompact;
                        brent_max_cellC_ = curCellInfo.centre();
                        brent_max_cellAlpha_ = pAlpha1[faceI];
                        brent_max_cell_nHat_ = pnHat[faceI];
                        brent_max_cell_gradAlpha1_ = pGradAlpha1[faceI];
                    }
                    
                    if(intfc_face.size())
                    {
                        pC_intfc[faceI] = intfc_face.centre(pts);
                        pA_intfc[faceI] = intfc_face.mag(pts);
                    }
                    else
                    {
                        pC_intfc[faceI] = curCellInfo.centre();
                        pA_intfc[faceI] = 0;
                    }

                    pC_ph1[faceI] = ph1_cell.centre(pts, fcs);
                    pC_ph0[faceI] = ph0_cell.centre(pts, fcs);
             
                    if(debugIR_)
                    {
                        os<< "Phase 1 sub-cell:" << endl;
                        Foam::plicFuncs::display_cell(ph1_cell, fcs, pts);
                        os<< "Phase 0 sub-cell:" << endl;
                        Foam::plicFuncs::display_cell(ph0_cell, fcs, pts);    
                        os<< "Interface face:" << endl;
                        Foam::plicFuncs::display_face(intfc_face, pts);
                        os<< "Interface centre:      " << pC_intfc[faceI] << nl
                            << "Interface area:        " << pA_intfc[faceI] << nl
                            << "Phase 1 centroid:      " << pC_ph1[faceI] << nl
                            << "Phase 0 centroid:      " << pC_ph0[faceI] << nl
                            << "-------------------------------------------------------" << nl
                            << "        Phase values on boundary face " << nfCompact <<  nl
                            << "-------------------------------------------------------" << nl 
                            << endl;
                    }
   
                    if(ph1_fcLbls[0] == -1)
                    {
                        face ph0_curFc = fcs[ph0_fcLbls[0]];
                        Cf_ph0_nei_[nfCompact] = ph0_curFc.centre(pts);
                        Af_ph0_nei_[nfCompact] = ph0_curFc.mag(pts);
                        Cf_ph1_nei_[nfCompact] = Cf_ph0_nei_[nfCompact];
                        Af_ph1_nei_[nfCompact] = 0;
                        face_phaseState_nei_[nfCompact] = 0;
                    }//end if(ph1_fcLbls[0] == -1)
                    else if(ph0_fcLbls[0] == -1)
                    {
                        face ph1_curFc = fcs[ph1_fcLbls[0]];
                        Cf_ph1_nei_[nfCompact] = ph1_curFc.centre(pts);
                        Af_ph1_nei_[nfCompact] = ph1_curFc.mag(pts);
                        Cf_ph0_nei_[nfCompact] = Cf_ph1_nei_[nfCompact];
                        Af_ph0_nei_[nfCompact] = 0;
                        face_phaseState_nei_[nfCompact] = 1;
                    }//end if(ph1_fcLbls[0] == -1)
                    else if(ph1_fcLbls[0] == -2)
                    {
                        Cf_ph1_nei_[nfCompact] = faces[nfCompact].centre(points);
                        Af_ph1_nei_[nfCompact] = 0;
                        Cf_ph0_nei_[nfCompact] = faces[nfCompact].centre(points);
                        Af_ph0_nei_[nfCompact] = 0;
                        face_phaseState_nei_[nfCompact] = 3;
                    }//end if(ph1_fcLbls[0] == -1)
                    else
                    {
                        face ph0_curFc = fcs[ph0_fcLbls[0]];
                        face ph1_curFc = fcs[ph1_fcLbls[0]];                        
                        Cf_ph0_nei_[nfCompact] = ph0_curFc.centre(pts);
                        Af_ph0_nei_[nfCompact] = ph0_curFc.mag(pts);
                        Cf_ph1_nei_[nfCompact] = ph1_curFc.centre(pts);
                        Af_ph1_nei_[nfCompact] = ph1_curFc.mag(pts);
                        face_phaseState_nei_[nfCompact] = 2;
                    }//end if(ph1_fcLbls[0] == -1)                    
                }//end if((pAlpha1[faceI] < ALPHA_2PH_MAX) && (pAlpha1[faceI] > ALPHA_2PH_MIN) && (mag(pGradAlpha1[faceI]) > GRADALPHA_MIN))
                else
                {
                    pC_ph1[faceI] = curCellInfo.centre();
                    pC_ph0[faceI] = curCellInfo.centre();
                
                    Cf_ph0_nei_[nfCompact] = pMeshCf[faceI];
                    Cf_ph1_nei_[nfCompact] = pMeshCf[faceI];
                    
                    if(pAlpha1[faceI] <= ALPHA_2PH_MIN)
                    {
                        Af_ph0_nei_[nfCompact] = pMeshMagSf[faceI];
                        Af_ph1_nei_[nfCompact] = 0;
                        face_phaseState_nei_[nfCompact] = 0;
                    }
                    else
                    {
                        Af_ph0_nei_[nfCompact] = 0;
                        Af_ph1_nei_[nfCompact] = pMeshMagSf[faceI];
                        face_phaseState_nei_[nfCompact] = 1;
                    }                    
                }//end if((pAlpha1[faceI] < ALPHA_2PH_MAX) && (pAlpha1[faceI] > ALPHA_2PH_MIN) && (mag(pGradAlpha1[faceI]) > GRADALPHA_MIN))

                if(debugIR_)
                {                    
                    os<< "Phase 1 centroid: " << Cf_ph1_nei_[nfCompact] << "    Phase 0 centroid: " << Cf_ph0_nei_[nfCompact] << nl    
                        << "Phase 1 area: " << Af_ph1_nei_[nfCompact] << "    Phase 0 area: " << Af_ph0_nei_[nfCompact] << nl
                        << "Face phase state:  " << face_phaseState_nei_[nfCompact] << nl
                        << "-------------------------------------------------------" << nl 
                        << endl;
                }

                if(debugIR_)
                {
                    os<< "=======================================================" << endl;
                    os<< " Done Face index: " << nfCompact << "    Flat fld face index: " << nCompact << endl;
                    os<< "=======================================================" << endl << endl;
                }

                nfCompact++;
                nCompact++;
            }//end forAll(pAlpha1, faceI)            
        }//end if(isA<zeroGradientFvPatchScalarField>(pAlpha1))
        else if(pAlpha1.coupled())
        {
            label bndFaceI = nfCompact - mesh().nInternalFaces();

            forAll(pAlpha1, faceI)
            {
                Cf_ph1_bndNei[bndFaceI] = Cf_ph1_own_[nfCompact];
                Cf_ph0_bndNei[bndFaceI] = Cf_ph0_own_[nfCompact];
                Af_ph1_bndNei[bndFaceI] = Af_ph1_own_[nfCompact];
                Af_ph0_bndNei[bndFaceI] = Af_ph0_own_[nfCompact];
                face_phaseState_bndNei[bndFaceI] = face_phaseState_own_[nfCompact];
                nfCompact++;
                bndFaceI++;
            }
        }
        else if(pAlpha1.fixesValue())
        {
            forAll(pAlpha1, faceI)
            {
                const cellInfo& curCellInfo = cells_flatFld_[nCompact];
                cell curCell(6);
                for(label fcI=0; fcI<curCell.size(); fcI++)
                {
                    curCell[fcI] = fcI;
                } 

                if(debugIR_)
                {
                    os<< "=======================================================" << nl
                        << "  Face index: " << nfCompact << "    Flat fld face index: " << nCompact << nl
                        << "=======================================================" << nl
                        << "Boundary alpha1 value:  " << pAlpha1[faceI] << nl                        
                        << "Patch face boundary cell: " << endl;
                    Foam::plicFuncs::display_cell(curCell, curCellInfo.faces(), curCellInfo.points());                        
                }

                pC_ph1[faceI] = curCellInfo.centre();
                pC_ph0[faceI] = curCellInfo.centre();
                
                Cf_ph0_nei_[nfCompact] = pMeshCf[faceI];
                Cf_ph1_nei_[nfCompact] = pMeshCf[faceI];
                    
                if(pAlpha1[faceI] <= ALPHA_2PH_MIN)
                {
                    Af_ph0_nei_[nfCompact] = pMeshMagSf[faceI];
                    Af_ph1_nei_[nfCompact] = 0;
                    face_phaseState_nei_[nfCompact] = 0;
                }
                else
                {
                    Af_ph0_nei_[nfCompact] = 0;
                    Af_ph1_nei_[nfCompact] = pMeshMagSf[faceI];
                    face_phaseState_nei_[nfCompact] = 1;
                }                    

                if(debugIR_)
                {                    
                    os<< "Phase 1 centroid: " << Cf_ph1_nei_[nfCompact] << "    Phase 0 centroid: " << Cf_ph0_nei_[nfCompact] << nl    
                        << "Phase 1 area: " << Af_ph1_nei_[nfCompact] << "    Phase 0 area: " << Af_ph0_nei_[nfCompact] << nl
                        << "Face phase state:  " << face_phaseState_nei_[nfCompact] << nl
                        << "-------------------------------------------------------" << nl 
                        << endl;
                }

                if(debugIR_)
                {
                    os<< "=======================================================" << endl;
                    os<< " Done Face index: " << nfCompact << "    Flat fld face index: " << nCompact << endl;
                    os<< "=======================================================" << endl << endl;
                }

                nfCompact++;
                nCompact++;
            }//end forAll(pAlpha1, faceI)
        }//end else if(isA<fixedValueFvPatchScalarField>(pAlpha1))
        else
        {
            forAll(pAlpha1, faceI)
            {
                const cellInfo& curCellInfo = cells_flatFld_[nCompact];
                cell curCell(6);
                for(label fcI=0; fcI<curCell.size(); fcI++)
                {
                    curCell[fcI] = fcI;
                } 

                if(debugIR_)
                {
                    os<< "=======================================================" << nl
                        << "  Face index: " << nfCompact << "    Flat fld face index: " << nCompact << nl
                        << "=======================================================" << nl
                        << "Boundary alpha1 value:  " << pAlpha1[faceI] << nl                        
                        << "Patch face boundary cell: " << endl;
                    Foam::plicFuncs::display_cell(curCell, curCellInfo.faces(), curCellInfo.points());                        
                }

                pC_ph1[faceI] = curCellInfo.centre();
                pC_ph0[faceI] = curCellInfo.centre();

                nfCompact++;
                nCompact++;
            }//end forAll(pAlpha1, faceI)
        }//end if(isA<zeroGradientFvPatchScalarField>(pAlpha1))

        if(debugIR_)
        {
            os<< "=======================================================" << endl;
            os<< "    Done Patch name:  " << patchNames[patchI] << endl;
            os<< "=======================================================" << endl << endl;                
        }
    }//end forAll(alpha_ph1_.boundaryField(), patchI)   

    syncTools::swapBoundaryFaceList(mesh(), Cf_ph1_bndNei);
    syncTools::swapBoundaryFaceList(mesh(), Cf_ph0_bndNei);
    syncTools::swapBoundaryFaceList(mesh(), Af_ph1_bndNei);
    syncTools::swapBoundaryFaceList(mesh(), Af_ph0_bndNei);
    syncTools::swapBoundaryFaceList(mesh(), face_phaseState_bndNei);

    forAll(alpha_ph1_.boundaryField(), patchI)
    {
        const polyPatch& pp = patches[patchI];

        label nfCompact = pp.start();
        label bndFaceI = nfCompact - mesh().nInternalFaces();

        if(pp.coupled())
        {
            forAll(pp, faceI)
            {
                Cf_ph1_nei_[nfCompact] = Cf_ph1_bndNei[bndFaceI];
                Cf_ph0_nei_[nfCompact] = Cf_ph0_bndNei[bndFaceI];
                Af_ph1_nei_[nfCompact] = Af_ph1_bndNei[bndFaceI];
                Af_ph0_nei_[nfCompact] = Af_ph0_bndNei[bndFaceI];
                face_phaseState_nei_[nfCompact] = face_phaseState_bndNei[bndFaceI];
                nfCompact++;
                bndFaceI++;
            }
        }
    }

    if(debugIR_)
    {
        os<< "//====================================================================\\" << nl
            << "          Done interface plane reconstruction in patch face cells" << nl
            << "\\====================================================================//" << nl
            << endl;
    }

    //os<< "Max Brent iterations: " << brent_iters_max << nl << endl;
    os<< mesh().time().timeName() << nl
        << "Max Brent iters = " << brent_iters_max_
        << "  Max Brent error = " << brent_err_max_ << endl;
        //<< "Max Brent iters cell: " << brent_max_cell_ << nl
        //<< "Max Brent iters cell loc: " << brent_max_cellC_ << nl
        //<< "Max Brent iters cell alpha1: " << brent_max_cellAlpha_ << nl
        //<< "Max Brent iters cell intfc normal: " << brent_max_cell_nHat_ << nl
        //<< "Max Brent iters cell grad alpha1: " << brent_max_cell_gradAlpha1_ << endl;
}


void Foam::plic::calc_face_phaseFluxes()
{    
    const polyBoundaryMesh& patches = mesh().boundaryMesh();

    if(debugF_)
    {
        Info<< "//====================================================================\\" << nl
            << "                 Interpolating U field to mesh points" << nl
            << "\\====================================================================//" << nl << endl;
    }

    //U_pts_ = ptInterp_.interpolate(U_);
    ptInterp_.interpolate(U_, U_pts_);

    if(debugF_)
    {
        Foam::plicFuncs::write_point_field(U_pts_, mesh());
    }

    if(debugF_)
    {
        Info<< "//====================================================================\\" << nl
            << "                 Done interpolating U field to mesh points" << nl
            << "\\====================================================================//" << nl << endl;
    }

    if(debugF_)
    {
        Info<< "//====================================================================\\" << nl
            << "               Collecting and distributing interface info" << nl
            << "\\====================================================================//" << nl << endl;
    }

    intfcInfo_collectData();

    if(debugF_)
    {
        Info<< "//====================================================================\\" << nl
            << "            Done collecting and distributing interface info" << nl
            << "\\====================================================================//" << nl << endl;
    }

    if(debugF2_)
    {
        Info<< "//====================================================================\\" << nl
            << "                  Tagging near-interface cells" << nl
            << "\\====================================================================//" << nl << nl
            << "------------------------------------------------------------------------" << nl
            << "cell index        mag(grad(alpha1))            near-intfc?" << nl
            << "------------------------------------------------------------------------" << endl;
    }

    forAll(cell_near_intfc_, cellI)
    {
        if(mag(gradAlpha1_.internalField()[cellI]) > GRADALPHA_MIN)
        {
            cell_near_intfc_[cellI] = true;
        }
        else
        {
            cell_near_intfc_[cellI] = false;
        }

        if(debugF2_)
        {
            Info<< "    " << cellI << "                  " << mag(gradAlpha1_.internalField()[cellI]) << "                              " << cell_near_intfc_[cellI] << endl; 
        }
    }
    
    if(debugF_)
    {
        Info<< nl
            << "//====================================================================\\" << nl
            << "                Done tagging near-interface cells" << nl
            << "\\====================================================================//" << nl << endl;
    }

    if(debugF2_)
    {
        Info<< "//====================================================================\\" << nl
            << "                 Tagging faces requiring 2-phase flux" << nl
            << "\\====================================================================//" << nl << nl
            << "========================================================================" << nl
            << "                           Internal faces" << nl
            << "========================================================================" << nl << nl
            << "------------------------------------------------------------------------" << nl
            << "  face      own    own ni?      nei    nei ni?      need 2-phase flux?" << nl
            << "------------------------------------------------------------------------" << endl;
    }

    const labelList& own = mesh().faceOwner();
    const labelList& nei = mesh().faceNeighbour();

    for(label faceI=0; faceI<mesh().nInternalFaces(); faceI++)
    {
        if(cell_near_intfc_[own[faceI]] || cell_near_intfc_[nei[faceI]])
        {
            fc_2ph_flux_needed_[faceI] = true;
        }
        else
        {
            fc_2ph_flux_needed_[faceI] = false;
        }
        
        if(debugF2_)
        {
            Info<< "  " << faceI << "         " << own[faceI] << "      " << cell_near_intfc_[own[faceI]] << "          " << nei[faceI] << "      " << cell_near_intfc_[nei[faceI]]
                << "                  " << fc_2ph_flux_needed_[faceI] << endl;
        }
    }

    if(debugF2_)
    {
        Info<< nl
            << "========================================================================" << nl
            << "                           Boundary faces" << nl
            << "========================================================================" << nl << endl;            
    }

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        const fvPatchVectorField& pGradAlpha1 = gradAlpha1_.boundaryField()[patchI];
        label faceI = pGradAlpha1.patch().start();

        if(debugF2_)
        {
            Info<< "--------------------------------------------------------------------" << nl
                << "                         patch: " << patchI << nl
                << "--------------------------------------------------------------------" << nl << endl;
        }

        if(pp.coupled())
        {
            if(debugF2_)
            {
                Info<< "Coupled patch..." << nl << nl
                    << "------------------------------------------------------------------------" << nl
                    << "  face          ownGrad          neiGrad          need 2-phase flux?" << nl
                    << "------------------------------------------------------------------------" << endl;
            }

            vectorField nHatOwn(pGradAlpha1.patchInternalField());
            vectorField nHatNei(pGradAlpha1.patchNeighbourField());
            forAll(pGradAlpha1, i)
            {
                scalar gradAlphaOwn = mag(nHatOwn[i]);
                scalar gradAlphaNei = mag(nHatNei[i]);

                if(gradAlphaOwn > GRADALPHA_MIN || gradAlphaNei > GRADALPHA_MIN)
                {
                    fc_2ph_flux_needed_[faceI] = true;
                }
                else
                {
                    fc_2ph_flux_needed_[faceI] = false;
                }

                if(debugF2_)
                {
                    Info<< "  " << faceI << "        " << gradAlphaOwn << "        " << gradAlphaNei << "        " << fc_2ph_flux_needed_[faceI] << endl;
                }

                faceI++;
            }
        }
        else if(isA<emptyPolyPatch>(pp))
        {
            if(debugF2_)
            {
                Info<< "Empty patch..." << nl << endl;
            }

            forAll(pGradAlpha1, i)
            {
                fc_2ph_flux_needed_[faceI] = false;
                faceI++;
            }
        }
        else
        {
            if(debugF2_)
            {
                Info<< "Non-coupled patch..." << nl << nl
                    << "------------------------------------------------------------------------" << nl
                    << "  face        own      own ni?          need 2-phase flux?" << nl
                    << "------------------------------------------------------------------------" << endl;
            }

            forAll(pGradAlpha1, i)
            {                
                if(cell_near_intfc_[own[faceI]])
                {
                    fc_2ph_flux_needed_[faceI] = true;
                }
                else
                {
                    fc_2ph_flux_needed_[faceI] = false;
                }

                if(debugF2_)
                {
                    Info<< "  " << faceI << "           " << own[faceI] << "        " << cell_near_intfc_[own[faceI]] << "                " << fc_2ph_flux_needed_[faceI] << endl;
                }

                faceI++;
            }
        }
    }

    if(debugF_)
    {
        Info<< "//====================================================================\\" << nl
            << "                 Done tagging faces requiring 2-phase flux" << nl
            << "\\====================================================================//" << nl << nl << endl;
    }

    if(debugF_)
    {
        Info<< "//====================================================================\\" << nl
            << "                  Calculating face 2-phase fluxes" << nl
            << "\\====================================================================//" << nl << nl << endl;
    }

    List<tetPoints> curFaceFluxTets(20);
    List<scalar> curFaceFluxTetVols(20);

    surfaceScalarField fAlpha1(fvc::interpolate(alpha_ph1_,"alpha1"));

    if(debugF_)
    {
        Info<< "========================================================================" << nl
            << "                    Internal face 2-phase fluxes" << nl
            << "========================================================================" << nl << endl;
    }

    for(label faceI=0; faceI<mesh().nInternalFaces(); faceI++)
    {
        if(fc_2ph_flux_needed_[faceI] && (mag(phi_[faceI]) > SMALLEST_PHI))
        {
            if(debugF_)
            {
                Info<< "========================================================================" << nl
                    << "               Calculating 2-phase flux for face " << faceI << nl
                    << "========================================================================" << nl << endl;
            }

            scalar curFaceFlux = 0;
            
            if(debugF_)
            {
                Info<< "========================================================================" << nl
                    << "              Step 1: Calculating face flux tetrahedrons " << nl
                    << "========================================================================" << nl << endl;
            }

            scalar curFacePhi = phi_[faceI];

            calcFaceFluxTets(faceI, curFacePhi, curFaceFluxTets, curFaceFluxTetVols);

            if(debugF_)
            {
                Info<< "========================================================================" << nl
                    << "            Done Step 1: Calculating face flux tetrahedrons" << nl
                    << "========================================================================" << nl << endl;
            }
            
            // clip each of 20 tets with the cells in the stencil for the
            // current face to find from which cell that portion of flux
            // comes from                
            const labelList& curFaceCells = faceStencil().stencil()[faceI];

            if(debugF_)
            {
                Info<< "Stencil for face " << faceI << endl;
                Foam::plicFuncs::display_labelList(curFaceCells);
            }

            if(debugF_)
            {
                Info<< "========================================================================" << nl
                    << "   Step 2: Intersecting face flux tetrahedrons with phase-1 subcells" << nl
                    << "========================================================================" << nl << endl;
            }

            for(label tetI=0; tetI<curFaceFluxTets.size(); tetI++)
            {
                if(debugF_)
                {
                    Info<< "--------------------------------------------------------------------" << nl
                        << "           Calculating contributions of tetrahedron " << tetI << nl
                        << "--------------------------------------------------------------------" << nl << endl;
                }

                scalar curTetVolFlux = 0;
                for(label cellI=0; cellI<curFaceCells.size(); cellI++)
                {
                    if(debugF_)
                    {
                        Info<< "--------------------------------------------------------------------" << nl
                            << "                   Intersecting with cell " << curFaceCells[cellI] << nl
                            << "--------------------------------------------------------------------" << nl << endl;
                    }
                    curTetVolFlux += tet_cell_intersect(curFaceFluxTets[tetI], curFaceCells[cellI]);
                }
                curFaceFlux += curTetVolFlux*sign(curFaceFluxTetVols[tetI]);

                if(debugF_)
                {
                    Info<< "Terahedron " << tetI << nl
                        << "Flux contributed:  " << curTetVolFlux << nl
                        << "Sign of flux:  " << sign(curFaceFluxTetVols[tetI]) << nl
                        << "Total face flux:  " << curFaceFlux << nl << endl;
                }
            }

            if(debugF_)
            {
                Info<< "========================================================================" << nl
                    << " Done Step 2: Intersecting face flux tetrahedrons with phase-1 subcells" << nl
                    << "========================================================================" << nl << endl;
            }

            phiAlpha1_[faceI] = curFaceFlux/mesh().time().deltaT().value();

            if(debugF_)
            {
                Info<< "Calculated face flux of phase-1:  " << phiAlpha1_[faceI] << nl << endl;
            }

            if(debugF_)
            {
                Info<< "========================================================================" << nl
                    << "              Done calculating 2-phase flux for face " << faceI << nl
                    << "========================================================================" << nl << endl;
            }
        }
        else
        {
            if(debugF_)
            {
                Info<< "========================================================================" << nl
                    << "               Calculating 1-phase flux for face " << faceI << nl
                    << "========================================================================" << nl << endl;
            }

            phiAlpha1_[faceI] = phi_[faceI]*fAlpha1[faceI];

            if(debugF_)
            {
                Info<< "Calculated face flux of phase-1:  " << phiAlpha1_[faceI] << nl << endl;
            }

            if(debugF_)
            {
                Info<< "========================================================================" << nl
                    << "             Done calculating 1-phase flux for face " << faceI << nl
                    << "========================================================================" << nl << endl;
            }
        }
    }

    if(debugF_)
    {
        Info<< nl
            << "========================================================================" << nl
            << "                    Boundary face 2-phase fluxes" << nl
            << "========================================================================" << nl << endl;
    }

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        const fvsPatchScalarField& pphi = phi_.boundaryField()[patchI];
        fvsPatchScalarField& pphiAlpha1 = phiAlpha1_.boundaryField()[patchI];
        fvsPatchScalarField& pfAlpha1 = fAlpha1.boundaryField()[patchI];
        label faceI = pp.start();

        if(debugF_)
        {
            Info<< "--------------------------------------------------------------------" << nl
                << "                         patch: " << patchI << nl
                << "--------------------------------------------------------------------" << nl << endl;
        }

        forAll(pphi, i)
        {
            if(fc_2ph_flux_needed_[faceI] && (mag(pphi[i]) > SMALLEST_PHI))
            {
                if(debugF_)
                {
                    Info<< "========================================================================" << nl
                        << "               Calculating 2-phase flux for face " << faceI << nl
                        << "========================================================================" << nl << endl;
                }

                scalar curFaceFlux = 0;

                if(debugF_)
                {
                    Info<< "========================================================================" << nl
                        << "              Step 1: Calculating face flux tetrahedrons " << nl
                        << "========================================================================" << nl << endl;
                }
            
                scalar curFacePhi = pphi[i];

                calcFaceFluxTets(faceI, curFacePhi, curFaceFluxTets, curFaceFluxTetVols);
            
                if(debugF_)
                {
                    Info<< "========================================================================" << nl
                        << "            Done Step 1: Calculating face flux tetrahedrons" << nl
                        << "========================================================================" << nl << endl;
                }

                // clip each of 20 tets with the cells in the stencil for the
                // current face to find from which cell that portion of flux
                // comes from    
                const labelList& curFaceCells = faceStencil().stencil()[faceI];

                if(debugF_)
                {
                    Info<< "Stencil for face " << faceI << endl;
                    Foam::plicFuncs::display_labelList(curFaceCells);
                }

                if(debugF_)
                {
                    Info<< "========================================================================" << nl
                        << "   Step 2: Intersecting face flux tetrahedrons with phase-1 subcells" << nl
                        << "========================================================================" << nl << endl;
                }

                for(label tetI=0; tetI<curFaceFluxTets.size(); tetI++)
                {
                    if(debugF_)
                    {
                        Info<< "--------------------------------------------------------------------" << nl
                            << "           Calculating contributions of tetrahedron " << tetI << nl
                            << "--------------------------------------------------------------------" << nl << endl;
                    }

                    scalar curTetVolFlux = 0;
                    for(label cellI=0; cellI<curFaceCells.size(); cellI++)
                    {
                        if(debugF_)
                        {
                            Info<< "--------------------------------------------------------------------" << nl
                                << "                   Intersecting with cell " << curFaceCells[cellI] << nl
                                << "--------------------------------------------------------------------" << nl << endl;
                        }

                        curTetVolFlux += tet_cell_intersect(curFaceFluxTets[tetI], curFaceCells[cellI]);
                    }
                    curFaceFlux += curTetVolFlux*sign(curFaceFluxTetVols[tetI]);

                    if(debugF_)
                    {
                        Info<< "Terahedron " << tetI << nl
                            << "Flux contributed:  " << curTetVolFlux << nl
                            << "Sign of flux:  " << sign(curFaceFluxTetVols[tetI]) << nl
                            << "Total face flux:  " << curFaceFlux << nl << endl;
                    }
                }

                if(debugF_)
                {
                    Info<< "========================================================================" << nl
                        << " Done Step 2: Intersecting face flux tetrahedrons with phase-1 subcells" << nl
                        << "========================================================================" << nl << endl;
                }

                pphiAlpha1[i] = curFaceFlux/mesh().time().deltaT().value();

                if(debugF_)
                {
                    Info<< "Calculated face flux of phase-1:  " << phiAlpha1_[faceI] << nl << endl;
                }

                if(debugF_)
                {
                    Info<< "========================================================================" << nl
                        << "              Done calculating 2-phase flux for face " << faceI << nl
                        << "========================================================================" << nl << endl;
                }
            }
            else
            {             
                if(debugF_)
                {
                    Info<< "========================================================================" << nl
                        << "               Calculating 1-phase flux for face " << faceI << nl
                        << "========================================================================" << nl << endl;
                }

                pphiAlpha1[i] = pphi[i]*pfAlpha1[i];

                if(debugF_)
                {
                    Info<< "Calculated face flux of phase-1:  " << phiAlpha1_[faceI] << nl << endl;
                }

                if(debugF_)
                {
                    Info<< "========================================================================" << nl
                        << "             Done calculating 1-phase flux for face " << faceI << nl
                        << "========================================================================" << nl << endl;
                }
            }
            faceI++;
        }
    }

    if(debugF_)
    {
        Info<< "//====================================================================\\" << nl
            << "                 Done calculating face 2-phase fluxes" << nl
            << "\\====================================================================//" << nl << nl << endl;
    }
}


void Foam::plic::calc_2ph_advFluxes
(
    const PtrList<volScalarField>& Y1,
    const PtrList<volScalarField>& Y0,
    PtrList<surfaceScalarField>& advFlux_Y1,
    PtrList<surfaceScalarField>& advFlux_Y0
)
{    
    const polyBoundaryMesh& patches = mesh().boundaryMesh();

    if(debugF_)
    {
        Info<< "//====================================================================\\" << nl
            << "                 Interpolating U field to mesh points" << nl
            << "\\====================================================================//" << nl << endl;
    }

    //U_pts_ = ptInterp_.interpolate(U_);
    ptInterp_.interpolate(U_, U_pts_);

    if(debugF2_)
    {
        Foam::plicFuncs::write_point_field(U_pts_, mesh());
    }

    if(debugF_)
    {
        Info<< "//====================================================================\\" << nl
            << "                 Done interpolating U field to mesh points" << nl
            << "\\====================================================================//" << nl << endl;
    }

    if(debugF_)
    {
        Info<< "//====================================================================\\" << nl
            << "               Collecting and distributing interface info" << nl
            << "\\====================================================================//" << nl << endl;
    }

    intfcInfo_collectData();

    Y_collectData(Y1, Y0);

    if(debugF_)
    {
        Info<< "//====================================================================\\" << nl
            << "            Done collecting and distributing interface info" << nl
            << "\\====================================================================//" << nl << endl;
    }

    if(debugF2_)
    {
        Info<< "//====================================================================\\" << nl
            << "                  Tagging near-interface cells" << nl
            << "\\====================================================================//" << nl << nl
            << "------------------------------------------------------------------------" << nl
            << "cell index        mag(grad(alpha1))            near-intfc?" << nl
            << "------------------------------------------------------------------------" << endl;
    }

    forAll(cell_near_intfc_, cellI)
    {
        if(mag(gradAlpha1_.internalField()[cellI]) > GRADALPHA_MIN)
        {
            cell_near_intfc_[cellI] = true;
        }
        else
        {
            cell_near_intfc_[cellI] = false;
        }

        if(debugF2_)
        {
            Info<< "    " << cellI << "                  " << mag(gradAlpha1_.internalField()[cellI]) << "                              " << cell_near_intfc_[cellI] << endl; 
        }
    }
    
    if(debugF_)
    {
        Info<< nl
            << "//====================================================================\\" << nl
            << "                Done tagging near-interface cells" << nl
            << "\\====================================================================//" << nl << endl;
    }

    if(debugF2_)
    {
        Info<< "//====================================================================\\" << nl
            << "                 Tagging faces requiring 2-phase flux" << nl
            << "\\====================================================================//" << nl << nl
            << "========================================================================" << nl
            << "                           Internal faces" << nl
            << "========================================================================" << nl << nl
            << "------------------------------------------------------------------------" << nl
            << "  face      own    own ni?      nei    nei ni?      need 2-phase flux?" << nl
            << "------------------------------------------------------------------------" << endl;
    }

    const labelList& own = mesh().faceOwner();
    const labelList& nei = mesh().faceNeighbour();

    for(label faceI=0; faceI<mesh().nInternalFaces(); faceI++)
    {
        if(cell_near_intfc_[own[faceI]] || cell_near_intfc_[nei[faceI]])
        {
            fc_2ph_flux_needed_[faceI] = true;
        }
        else
        {
            fc_2ph_flux_needed_[faceI] = false;
        }
        
        if(debugF2_)
        {
            Info<< "  " << faceI << "         " << own[faceI] << "      " << cell_near_intfc_[own[faceI]] << "          " << nei[faceI] << "      " << cell_near_intfc_[nei[faceI]]
                << "                  " << fc_2ph_flux_needed_[faceI] << endl;
        }
    }

    if(debugF2_)
    {
        Info<< nl
            << "========================================================================" << nl
            << "                           Boundary faces" << nl
            << "========================================================================" << nl << endl;            
    }

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        const fvPatchVectorField& pGradAlpha1 = gradAlpha1_.boundaryField()[patchI];
        label faceI = pGradAlpha1.patch().start();

        if(debugF2_)
        {
            Info<< "--------------------------------------------------------------------" << nl
                << "                         patch: " << patchI << nl
                << "--------------------------------------------------------------------" << nl << endl;
        }

        if(pp.coupled())
        {
            if(debugF2_)
            {
                Info<< "Coupled patch..." << nl << nl
                    << "------------------------------------------------------------------------" << nl
                    << "  face          ownGrad          neiGrad          need 2-phase flux?" << nl
                    << "------------------------------------------------------------------------" << endl;
            }

            vectorField nHatOwn(pGradAlpha1.patchInternalField());
            vectorField nHatNei(pGradAlpha1.patchNeighbourField());
            forAll(pGradAlpha1, i)
            {
                scalar gradAlphaOwn = mag(nHatOwn[i]);
                scalar gradAlphaNei = mag(nHatNei[i]);

                if(gradAlphaOwn > GRADALPHA_MIN || gradAlphaNei > GRADALPHA_MIN)
                {
                    fc_2ph_flux_needed_[faceI] = true;
                }
                else
                {
                    fc_2ph_flux_needed_[faceI] = false;
                }

                if(debugF2_)
                {
                    Info<< "  " << faceI << "        " << gradAlphaOwn << "        " << gradAlphaNei << "        " << fc_2ph_flux_needed_[faceI] << endl;
                }

                faceI++;
            }
        }
        else if(isA<emptyPolyPatch>(pp))
        {
            if(debugF2_)
            {
                Info<< "Empty patch..." << nl << endl;
            }

            forAll(pGradAlpha1, i)
            {
                fc_2ph_flux_needed_[faceI] = false;
                faceI++;
            }
        }
        else
        {
            if(debugF2_)
            {
                Info<< "Non-coupled patch..." << nl << nl
                    << "------------------------------------------------------------------------" << nl
                    << "  face        own      own ni?          need 2-phase flux?" << nl
                    << "------------------------------------------------------------------------" << endl;
            }

            forAll(pGradAlpha1, i)
            {                
                if(cell_near_intfc_[own[faceI]])
                {
                    fc_2ph_flux_needed_[faceI] = true;
                }
                else
                {
                    fc_2ph_flux_needed_[faceI] = false;
                }

                if(debugF2_)
                {
                    Info<< "  " << faceI << "           " << own[faceI] << "        " << cell_near_intfc_[own[faceI]] << "                " << fc_2ph_flux_needed_[faceI] << endl;
                }

                faceI++;
            }
        }
    }

    if(debugF_)
    {
        Info<< "//====================================================================\\" << nl
            << "                 Done tagging faces requiring 2-phase flux" << nl
            << "\\====================================================================//" << nl << nl << endl;
    }

    if(debugF_)
    {
        Info<< "//====================================================================\\" << nl
            << "                  Calculating face 2-phase fluxes" << nl
            << "\\====================================================================//" << nl << nl << endl;
    }

    List<tetPoints> curFaceFluxTets(20);
    List<scalar> curFaceFluxTetVols(20);

    fAlpha1_ = fvc::interpolate(alpha_ph1_, "alpha1");
    frho1_ = fvc::interpolate(rho1_, "rho");
    frho0_ = fvc::interpolate(rho0_, "rho");

    for(label i=0; i<(nSpecies_ - 1); i++)
    {
        fY1_[i] = fvc::interpolate(Y1[i], "Yi");
        fY0_[i] = fvc::interpolate(Y0[i], "Yi");
        //plicFuncs::display_surfaceField(fY0_[i], mesh());

        /*
        forAll(fY0_[i].boundaryField(), patchI)
        {
            const polyPatch& pp = mesh().boundaryMesh()[patchI];
            fvsPatchScalarField& pfY0i = fY0_[i].boundaryField()[patchI];
            fvsPatchScalarField& pfY1i = fY0_[i].boundaryField()[patchI];
            const fvPatchScalarField& pY0i = Y0[i].boundaryField()[patchI];
            const fvPatchScalarField& pY1i = Y1[i].boundaryField()[patchI];
            if(isA<zeroGradientFvPatchScalarField>(pY0i))
            {
                label faceI = pp.start();

                forAll(pY0i, fcI)
                {
                    label faceOwn = own[faceI];
                    Info<< "face " << faceI
                        << nl
                        << "pfY0i[fcI]: " << pfY0i[fcI]
                        << "  Y0i[faceOwn]: " << Y0[i].internalField()[faceOwn]
                        << "  pY0i[fcI]: " << pY0i[fcI]
                        << endl;
                    pfY0i[fcI] = Y0[i].internalField()[faceOwn];
                    pfY1i[fcI] = Y1[i].internalField()[faceOwn];
                    faceI++;
                }
            }
        }
            */
    }

    if(debugF_)
    {
        Info<< "========================================================================" << nl
            << "                    Internal face 2-phase fluxes" << nl
            << "========================================================================" << nl << endl;
    }

    for(label faceI=0; faceI<mesh().nInternalFaces(); faceI++)
    {
        if(fc_2ph_flux_needed_[faceI] && (mag(phi_[faceI]) > SMALLEST_PHI))
        {
            if(debugF_)
            {
                Info<< "========================================================================" << nl
                    << "               Calculating 2-phase flux for face " << faceI << nl
                    << "========================================================================" << nl << endl;
            }

            scalar curFaceFlux_ph1 = 0;
            scalar curFaceFlux_ph0 = 0;
            List<scalar> curFaceY1Flux(nSpecies_ - 1);
            List<scalar> curFaceY0Flux(nSpecies_ - 1);
            for(label i=0; i<(nSpecies_ - 1); i++)
            {
                curFaceY1Flux[i] = 0;
                curFaceY0Flux[i] = 0;
            }
            
            if(debugF_)
            {
                Info<< "========================================================================" << nl
                    << "              Step 1: Calculating face flux tetrahedrons " << nl
                    << "========================================================================" << nl << endl;
            }

            scalar curFacePhi = phi_[faceI];

            calcFaceFluxTets(faceI, curFacePhi, curFaceFluxTets, curFaceFluxTetVols);

            if(debugF_)
            {
                Info<< "========================================================================" << nl
                    << "            Done Step 1: Calculating face flux tetrahedrons" << nl
                    << "========================================================================" << nl << endl;
            }
            
            // clip each of 20 tets with the cells in the stencil for the
            // current face to find from which cell that portion of flux
            // comes from                
            const labelList& curFaceCells = faceStencil().stencil()[faceI];

            if(debugF_)
            {
                Info<< "Stencil for face " << faceI << endl;
                Foam::plicFuncs::display_labelList(curFaceCells);
            }

            if(debugF_)
            {
                Info<< "========================================================================" << nl
                    << "   Step 2: Intersecting face flux tetrahedrons with phase-1 subcells" << nl
                    << "========================================================================" << nl << endl;
            }

            for(label tetI=0; tetI<curFaceFluxTets.size(); tetI++)
            {
                if(debugF_)
                {
                    Info<< "--------------------------------------------------------------------" << nl
                        << "           Calculating contributions of tetrahedron " << tetI << nl
                        << "--------------------------------------------------------------------" << nl << endl;
                }

                scalar curTetVolFlux_ph1 = 0;
                scalar curTetVolFlux_ph0 = 0;
                List<scalar> curTetY1Flux(nSpecies_ - 1);
                List<scalar> curTetY0Flux(nSpecies_ - 1);
                for(label i=0; i<(nSpecies_ - 1); i++)
                {
                    curTetY1Flux[i] = 0;
                    curTetY0Flux[i] = 0;
                }

                for(label cellI=0; cellI<curFaceCells.size(); cellI++)
                {
                    if(debugF_)
                    {
                        Info<< "--------------------------------------------------------------------" << nl
                            << "                   Intersecting with cell " << curFaceCells[cellI] << nl
                            << "--------------------------------------------------------------------" << nl << endl;
                    }

                    label curCell_lbl = curFaceCells[cellI];
                    scalar curTetCurCellVolFlux_ph1 = 0;
                    scalar curTetCurCellVolFlux_ph0 = 0;
                    tet_cell_intersect(curFaceFluxTets[tetI], curCell_lbl, curTetCurCellVolFlux_ph1, curTetCurCellVolFlux_ph0);
                    curTetVolFlux_ph1 += curTetCurCellVolFlux_ph1;
                    curTetVolFlux_ph0 += curTetCurCellVolFlux_ph0;

                    if(debugF_)
                    {
                        Info<< "curCell alpha1 = " << alpha_ph1_flatFld_[curCell_lbl] << nl
                            << "curTetCurCellVolFlux_ph1 = " << curTetCurCellVolFlux_ph1 << "  curTetCurCellVolFlux_ph0 = " << curTetCurCellVolFlux_ph0 << nl
                            << "curCell rho1 = " << rho1_flatFld_[curCell_lbl] << "  curCell rho0 = " << rho0_flatFld_[curCell_lbl] << endl;
                    }

                    for(label i=0; i<(nSpecies_ - 1); i++)
                    {
                        curTetY1Flux[i] += rho1_flatFld_[curCell_lbl]*Y1_flatFld_[i][curCell_lbl]*curTetCurCellVolFlux_ph1;
                        curTetY0Flux[i] += rho0_flatFld_[curCell_lbl]*Y0_flatFld_[i][curCell_lbl]*curTetCurCellVolFlux_ph0;
                        
                        if(debugF_)
                        {
                            Info<< "  Y1[" << i << "] = " << Y1_flatFld_[i][curCell_lbl] << "  curTet curCell Y1[" << i << "] flux = " << rho1_flatFld_[curCell_lbl]*Y1_flatFld_[i][curCell_lbl]*curTetCurCellVolFlux_ph1 << nl
                                << "  Y0[" << i << "] = " << Y0_flatFld_[i][curCell_lbl] << "  curTet curCell Y0[" << i << "] flux = " << rho0_flatFld_[curCell_lbl]*Y0_flatFld_[i][curCell_lbl]*curTetCurCellVolFlux_ph0 << endl;
                        }
                    }// end for(label i=0; i<(nSpecies_ - 1); i++)

                    if(debugF_)
                    {
                        Info<< nl
                            << "--------------------------------------------------------------------" << endl;                            
                    }
                }// end for(label cellI=0; cellI<curFaceCells.size(); cellI++)

                curFaceFlux_ph1 += curTetVolFlux_ph1*sign(curFaceFluxTetVols[tetI]);
                curFaceFlux_ph0 += curTetVolFlux_ph0*sign(curFaceFluxTetVols[tetI]);

                if(debugF_)
                {
                    Info<< "Tetrahedron " << tetI << nl
                        << "Phase-1 flux contribution:  " << curTetVolFlux_ph1 << nl
                        << "Phase-0 flux contribution:  " << curTetVolFlux_ph0 << nl
                        << "Phase-1 + Phase-0 flux:  " << curTetVolFlux_ph1 + curTetVolFlux_ph0 << nl
                        << "Tetrahedron vol:  " << curFaceFluxTetVols[tetI] << nl
                        << "Tetrahedron vol flux error: " << mag(curTetVolFlux_ph1 + curTetVolFlux_ph0 - mag(curFaceFluxTetVols[tetI])) << nl
                        << "Sign of flux:  " << sign(curFaceFluxTetVols[tetI]) << nl
                        << "Cumulative phase-1 face flux:  " << curFaceFlux_ph1 << nl
                        << "Cumulative phase-0 face flux:  " << curFaceFlux_ph0 << endl;

                    if(mag(curTetVolFlux_ph1 + curTetVolFlux_ph0 - mag(curFaceFluxTetVols[tetI])) > 1E-15)
                    {
                        Info<< "Tetrahedron vol flux error > 1E-15!!"
                            << endl;
                    }
                }

                for(label i=0; i<(nSpecies_ - 1); i++)
                {
                    curFaceY1Flux[i] += curTetY1Flux[i]*sign(curFaceFluxTetVols[tetI]);
                    curFaceY0Flux[i] += curTetY0Flux[i]*sign(curFaceFluxTetVols[tetI]);

                    if(debugF_)
                    {
                        Info<< "Y1[" << i << "] flux contribution:  " << curTetY1Flux[i] << nl                                              
                            << "Cumulative Y1[" << i << "] face flux:  " << curFaceY1Flux[i] << nl
                            << "Y0[" << i << "] flux contribution:  " << curTetY0Flux[i] << nl                                              
                            << "Cumulative Y0[" << i << "] face flux:  " << curFaceY0Flux[i] << endl;                            
                    }
                }// end for(label i=0; i<(nSpecies_ - 1); i++)

                if(debugF_)
                {
                    Info<< nl
                        << "--------------------------------------------------------------------" << endl;                            
                }
            }// end for(label tetI=0; tetI<curFaceFluxTets.size(); tetI++)

            if(debugF_)
            {
                Info<< "========================================================================" << nl
                    << " Done Step 2: Intersecting face flux tetrahedrons with phase-1 subcells" << nl
                    << "========================================================================" << nl << endl;
            }

            phiAlpha1_[faceI] = curFaceFlux_ph1/mesh().time().deltaT().value();
            phiAlpha0_[faceI] = curFaceFlux_ph0/mesh().time().deltaT().value();

            if(debugF_)
            {
                Info<< "Total phase-1 face flux:  " << phiAlpha1_[faceI] << nl
                    << "Total phase-0 face flux:  " << phiAlpha0_[faceI] << nl
                    << "Total phase-1 + phase-0 face flux:  " << phiAlpha1_[faceI] + phiAlpha0_[faceI] << nl
                    << "face phi:  " << curFacePhi << nl
                    << "face phi error:  " << mag(phiAlpha1_[faceI] + phiAlpha0_[faceI] - curFacePhi) << endl;

                if(mag(phiAlpha1_[faceI] + phiAlpha0_[faceI] - curFacePhi) > 1E-12)
                {
                    Info<< "face phi error > 1E-12!!"
                        << endl;
                }
            }

            for(label i=0; i<(nSpecies_ - 1); i++)
            {
                advFlux_Y1[i][faceI] = curFaceY1Flux[i]/mesh().time().deltaT().value();
                advFlux_Y0[i][faceI] = curFaceY0Flux[i]/mesh().time().deltaT().value();

                if(debugF_)
                {
                    Info<< "Total Y1[" << i << "] face flux:  " << advFlux_Y1[i][faceI] << nl
                        << "Total Y0[" << i << "] face flux:  " << advFlux_Y0[i][faceI] << endl;                        
                }
            }            

            if(debugF_)
            {
                Info<< nl
                    << "========================================================================" << nl
                    << "              Done calculating 2-phase flux for face " << faceI << nl
                    << "========================================================================" << nl << endl;
            }
        }// end if(fc_2ph_flux_needed_[faceI] && (mag(phi_[faceI]) > SMALLEST_PHI))
        else
        {
            if(debugF_)
            {
                Info<< "========================================================================" << nl
                    << "               Calculating 1-phase flux for face " << faceI << nl
                    << "========================================================================" << nl << endl;
            }

            phiAlpha1_[faceI] = phi_[faceI]*fAlpha1_[faceI];
            phiAlpha0_[faceI] = phi_[faceI]*(1 - fAlpha1_[faceI]);            

            if(debugF_)
            {
                Info<< "Face alpha1:  " << fAlpha1_[faceI] << nl
                    << "Total phase-1 face flux:  " << phiAlpha1_[faceI] << nl
                    << "Total phase-0 face flux:  " << phiAlpha0_[faceI] << nl
                    << "Total phase-1 + phase-0 face flux:  " << phiAlpha1_[faceI] + phiAlpha0_[faceI] << nl
                    << "face phi:  " << phi_[faceI] << nl
                    << "face phi error:  " << mag(phiAlpha1_[faceI] + phiAlpha0_[faceI] - phi_[faceI]) << endl;

                if(mag(phiAlpha1_[faceI] + phiAlpha0_[faceI] - phi_[faceI]) > 1E-12)
                {
                    Info<< "face phi error > 1E-12!!"
                        << endl;
                }
            }
            
            if(fAlpha1_[faceI] > 0.5)
            {
                for(label i=0; i<(nSpecies_ - 1); i++)
                {
                    advFlux_Y1[i][faceI] = phi_[faceI]*frho1_[faceI]*fY1_[i][faceI];
                    advFlux_Y0[i][faceI] = 0;

                    if(debugF_)
                    {
                        Info<< "Total Y1[" << i << "] face flux:  " << advFlux_Y1[i][faceI] << nl
                            << "Total Y0[" << i << "] face flux:  " << advFlux_Y0[i][faceI] << endl;                        
                    }
                }
            }// end if (fAlpha1_[faceI] > 0.5)
            else
            {
                for(label i=0; i<(nSpecies_ - 1); i++)
                {
                    advFlux_Y1[i][faceI] = 0;
                    advFlux_Y0[i][faceI] = phi_[faceI]*frho0_[faceI]*fY0_[i][faceI];

                    if(debugF_)
                    {
                        Info<< "Total Y1[" << i << "] face flux:  " << advFlux_Y1[i][faceI] << nl
                            << "Total Y0[" << i << "] face flux:  " << advFlux_Y0[i][faceI] << endl;                        
                    }
                }
            }// end if (fAlpha1_[faceI] > 0.5)

            if(debugF_)
            {
                Info<< "========================================================================" << nl
                    << "             Done calculating 1-phase flux for face " << faceI << nl
                    << "========================================================================" << nl << endl;
            }
        }// end !if(fc_2ph_flux_needed_[faceI] && (mag(phi_[faceI]) > SMALLEST_PHI))
    }

    if(debugF_)
    {
        Info<< nl
            << "========================================================================" << nl
            << "                    Boundary face 2-phase fluxes" << nl
            << "========================================================================" << nl << endl;
    }

    forAll(phi_.boundaryField(), patchI)
    {
        const polyPatch& pp = patches[patchI];
        const fvsPatchScalarField& pphi = phi_.boundaryField()[patchI];
        fvsPatchScalarField& pphiAlpha1 = phiAlpha1_.boundaryField()[patchI];
        fvsPatchScalarField& pphiAlpha0 = phiAlpha0_.boundaryField()[patchI];
        fvsPatchScalarField& pfAlpha1 = fAlpha1_.boundaryField()[patchI];
        label faceI = pp.start();

        if(debugF_)
        {
            Info<< "--------------------------------------------------------------------" << nl
                << "                         patch: " << patchI << nl
                << "--------------------------------------------------------------------" << nl << endl;
        }

        forAll(pphi, fcI)
        {
            if(fc_2ph_flux_needed_[faceI] && (mag(pphi[fcI]) > SMALLEST_PHI))
            {
                if(debugF_)
                {
                    Info<< "========================================================================" << nl
                        << "               Calculating 2-phase flux for face " << faceI << nl
                        << "========================================================================" << nl << endl;
                }
                
                scalar curFaceFlux_ph1 = 0;
                scalar curFaceFlux_ph0 = 0;                
                List<scalar> curFaceY1Flux((nSpecies_ - 1));
                List<scalar> curFaceY0Flux((nSpecies_ - 1));
                for(label i=0; i<(nSpecies_ - 1); i++)
                {
                    curFaceY1Flux[i] = 0;
                    curFaceY0Flux[i] = 0;
                }

                if(debugF_)
                {
                    Info<< "========================================================================" << nl
                        << "              Step 1: Calculating face flux tetrahedrons " << nl
                        << "========================================================================" << nl << endl;
                }
            
                scalar curFacePhi = pphi[fcI];

                calcFaceFluxTets(faceI, curFacePhi, curFaceFluxTets, curFaceFluxTetVols);
            
                if(debugF_)
                {
                    Info<< "========================================================================" << nl
                        << "            Done Step 1: Calculating face flux tetrahedrons" << nl
                        << "========================================================================" << nl << endl;
                }

                // clip each of 20 tets with the cells in the stencil for the
                // current face to find from which cell that portion of flux
                // comes from    
                const labelList& curFaceCells = faceStencil().stencil()[faceI];

                if(debugF_)
                {
                    Info<< "Stencil for face " << faceI << endl;
                    Foam::plicFuncs::display_labelList(curFaceCells);
                }

                if(debugF_)
                {
                    Info<< "========================================================================" << nl
                        << "   Step 2: Intersecting face flux tetrahedrons with phase-1 subcells" << nl
                        << "========================================================================" << nl << endl;
                }

                for(label tetI=0; tetI<curFaceFluxTets.size(); tetI++)
                {
                    if(debugF_)
                    {
                        Info<< "--------------------------------------------------------------------" << nl
                            << "           Calculating contributions of tetrahedron " << tetI << nl
                            << "--------------------------------------------------------------------" << nl << endl;
                    }

                    scalar curTetVolFlux_ph1 = 0;
                    scalar curTetVolFlux_ph0 = 0;
                    List<scalar> curTetY1Flux((nSpecies_ - 1));
                    List<scalar> curTetY0Flux((nSpecies_ - 1));
                    for(label i=0; i<(nSpecies_ - 1); i++)
                    {
                        curTetY1Flux[i] = 0;
                        curTetY0Flux[i] = 0;
                    }

                    for(label cellI=0; cellI<curFaceCells.size(); cellI++)
                    {
                        if(debugF_)
                        {
                            Info<< "--------------------------------------------------------------------" << nl
                                << "                   Intersecting with cell " << curFaceCells[cellI] << nl
                                << "--------------------------------------------------------------------" << nl << endl;
                        }
                        
                        label curCell_lbl = curFaceCells[cellI];
                        scalar curTetCurCellVolFlux_ph1 = 0;
                        scalar curTetCurCellVolFlux_ph0 = 0;
                        tet_cell_intersect(curFaceFluxTets[tetI], curCell_lbl, curTetCurCellVolFlux_ph1, curTetCurCellVolFlux_ph0);
                        curTetVolFlux_ph1 += curTetCurCellVolFlux_ph1;
                        curTetVolFlux_ph0 += curTetCurCellVolFlux_ph0;

                        if(debugF_)
                        {
                            Info<< "curCell alpha1 = " << alpha_ph1_flatFld_[curCell_lbl] << nl
                                << "curTetCurCellVolFlux_ph1 = " << curTetCurCellVolFlux_ph1 << "  curTetCurCellVolFlux_ph0 = " << curTetCurCellVolFlux_ph0 << nl
                                << "curCell rho1 = " << rho1_flatFld_[curCell_lbl] << "  curCell rho0 = " << rho0_flatFld_[curCell_lbl] << endl;
                        }

                        for(label i=0; i<(nSpecies_ - 1); i++)
                        {
                            curTetY1Flux[i] += rho1_flatFld_[curCell_lbl]*Y1_flatFld_[i][curCell_lbl]*curTetCurCellVolFlux_ph1;
                            curTetY0Flux[i] += rho0_flatFld_[curCell_lbl]*Y0_flatFld_[i][curCell_lbl]*curTetCurCellVolFlux_ph0;
                        
                            if(debugF_)
                            {
                                Info<< "  Y1[" << i << "] = " << Y1_flatFld_[i][curCell_lbl] << "  curTet curCell Y1[" << i << "] flux = " << rho1_flatFld_[curCell_lbl]*Y1_flatFld_[i][curCell_lbl]*curTetCurCellVolFlux_ph1 << nl
                                    << "  Y0[" << i << "] = " << Y0_flatFld_[i][curCell_lbl] << "  curTet curCell Y0[" << i << "] flux = " << rho0_flatFld_[curCell_lbl]*Y0_flatFld_[i][curCell_lbl]*curTetCurCellVolFlux_ph0 << endl;
                            }
                        }// end for(label i=0; i<(nSpecies_ - 1); i++)

                        if(debugF_)
                        {
                            Info<< nl
                                << "--------------------------------------------------------------------" << endl;                            
                        }
                    }// end for(label cellI=0; cellI<curFaceCells.size(); cellI++)

                    curFaceFlux_ph1 += curTetVolFlux_ph1*sign(curFaceFluxTetVols[tetI]);
                    curFaceFlux_ph0 += curTetVolFlux_ph0*sign(curFaceFluxTetVols[tetI]);                    

                    if(debugF_)
                    {
                        Info<< "Tetrahedron " << tetI << nl
                        << "Phase-1 flux contribution:  " << curTetVolFlux_ph1 << nl
                        << "Phase-0 flux contribution:  " << curTetVolFlux_ph0 << nl
                        << "Phase-1 + Phase-0 flux:  " << curTetVolFlux_ph1 + curTetVolFlux_ph0 << nl
                        << "Tetrahedron vol:  " << curFaceFluxTetVols[tetI] << nl
                        << "Tetrahedron vol flux error: " << mag(curTetVolFlux_ph1 + curTetVolFlux_ph0 - mag(curFaceFluxTetVols[tetI])) << nl
                        << "Sign of flux:  " << sign(curFaceFluxTetVols[tetI]) << nl
                        << "Cumulative phase-1 face flux:  " << curFaceFlux_ph1 << nl
                        << "Cumulative phase-0 face flux:  " << curFaceFlux_ph0 << endl;

                        if(mag(curTetVolFlux_ph1 + curTetVolFlux_ph0 - mag(curFaceFluxTetVols[tetI])) > 1E-15)
                        {
                            Info<< "Tetrahedron vol flux error > 1E-15!!"
                                << endl;
                        }
                    }

                    for(label i=0; i<(nSpecies_ - 1); i++)
                    {
                        curFaceY1Flux[i] += curTetY1Flux[i]*sign(curFaceFluxTetVols[tetI]);
                        curFaceY0Flux[i] += curTetY0Flux[i]*sign(curFaceFluxTetVols[tetI]);

                        if(debugF_)
                        {
                            Info<< "Y1[" << i << "] flux contribution:  " << curTetY1Flux[i] << nl                                              
                                << "Cumulative Y1[" << i << "] face flux:  " << curFaceY1Flux[i] << nl
                                << "Y0[" << i << "] flux contribution:  " << curTetY0Flux[i] << nl                                              
                                << "Cumulative Y0[" << i << "] face flux:  " << curFaceY0Flux[i] << endl;                            
                        }
                    }// end for(label i=0; i<(nSpecies_ - 1); i++)

                    if(debugF_)
                    {
                        Info<< nl
                            << "--------------------------------------------------------------------" << endl;                            
                    }
                }// end for(label tetI=0; tetI<curFaceFluxTets.size(); tetI++)

                if(debugF_)
                {
                    Info<< "========================================================================" << nl
                        << " Done Step 2: Intersecting face flux tetrahedrons with phase-1 subcells" << nl
                        << "========================================================================" << nl << endl;
                }

                pphiAlpha1[fcI] = curFaceFlux_ph1/mesh().time().deltaT().value();
                pphiAlpha0[fcI] = curFaceFlux_ph0/mesh().time().deltaT().value();

                if(debugF_)
                {
                    Info<< "Total phase-1 face flux:  " << pphiAlpha1[fcI] << nl
                        << "Total phase-0 face flux:  " << pphiAlpha0[fcI] << nl
                        << "phase-1 + phase-0 face flux:  " << pphiAlpha1[fcI] + pphiAlpha0[fcI] << nl
                        << "face phi:  " << curFacePhi << nl
                        << "face phi error:  " << mag(pphiAlpha1[fcI] + pphiAlpha0[fcI] - curFacePhi) << endl;

                    if(mag(pphiAlpha1[fcI] + pphiAlpha0[fcI] - curFacePhi) > 1E-12)
                    {
                        Info<< "face phi error > 1E-12!!"
                            << endl;
                    }
                }

                for(label i=0; i<(nSpecies_ - 1); i++)
                {
                    advFlux_Y1[i].boundaryField()[patchI][fcI] = curFaceY1Flux[i]/mesh().time().deltaT().value();
                    advFlux_Y0[i].boundaryField()[patchI][fcI] = curFaceY0Flux[i]/mesh().time().deltaT().value();

                    if(debugF_)
                    {
                        Info<< "Total Y1[" << i << "] face flux:  " << advFlux_Y1[i].boundaryField()[patchI][fcI] << nl
                            << "Total Y0[" << i << "] face flux:  " << advFlux_Y0[i].boundaryField()[patchI][fcI] << endl;                        
                    }
                }

                if(debugF_)
                {
                    Info<< "========================================================================" << nl
                        << "              Done calculating 2-phase flux for face " << faceI << nl
                        << "========================================================================" << nl << endl;
                }
            }// end if(fc_2ph_flux_needed_[faceI] && (mag(pphi[fcI]) > SMALLEST_PHI))
            else
            {             
                if(debugF_)
                {
                    Info<< "========================================================================" << nl
                        << "               Calculating 1-phase flux for face " << faceI << nl
                        << "========================================================================" << nl << endl;
                }

                pphiAlpha1[fcI] = pphi[fcI]*pfAlpha1[fcI];
                pphiAlpha0[fcI] = pphi[fcI]*(1 - pfAlpha1[fcI]);

                if(debugF_)
                {
                    Info<< "Face alpha1:  " << pfAlpha1[fcI] << nl
                        << "Total phase-1 face flux:  " << pphiAlpha1[fcI] << nl
                        << "Total phase-0 face flux:  " << pphiAlpha0[fcI] << nl
                        << "Total phase-1 + phase-0 face flux:  " << pphiAlpha1[fcI] + pphiAlpha0[fcI] << nl
                        << "face phi:  " << pphi[fcI] << nl
                        << "face phi error:  " << mag(pphiAlpha1[fcI] + pphiAlpha0[fcI] - pphi[fcI]) << endl;

                    if(mag(pphiAlpha1[fcI] + pphiAlpha0[fcI] - pphi[fcI]) > 1E-12)
                    {
                        Info<< "face phi error > 1E-12!!"
                            << endl;
                    }
                }

                if(pfAlpha1[fcI] > 0.5)
                {
                    for(label i=0; i<(nSpecies_ - 1); i++)
                    {
                        advFlux_Y1[i].boundaryField()[patchI][fcI] = pphi[fcI]*frho1_.boundaryField()[patchI][fcI]*fY1_[i].boundaryField()[patchI][fcI];
                        advFlux_Y0[i].boundaryField()[patchI][fcI] = 0;

                        if(debugF_)
                        {
                            Info<< "Total Y1[" << i << "] face flux:  " << advFlux_Y1[i].boundaryField()[patchI][fcI] << nl
                                << "Total Y0[" << i << "] face flux:  " << advFlux_Y0[i].boundaryField()[patchI][fcI] << endl;                        
                        }
                    }
                }// end if (pfAlpha1[fcI] > 0.5)
                else
                {
                    for(label i=0; i<(nSpecies_ - 1); i++)
                    {
                        advFlux_Y1[i].boundaryField()[patchI][fcI] = 0;
                        advFlux_Y0[i].boundaryField()[patchI][fcI] = pphi[fcI]*frho0_.boundaryField()[patchI][fcI]*fY0_[i].boundaryField()[patchI][fcI];

                        if(debugF_)
                        {
                            Info<< "Total Y1[" << i << "] face flux:  " << advFlux_Y1[i].boundaryField()[patchI][fcI] << nl
                                << "frho0: " << frho0_.boundaryField()[patchI][fcI] << "  fY0i: " <<  fY0_[i].boundaryField()[patchI][fcI] << nl
                                << "Total Y0[" << i << "] face flux:  " << advFlux_Y0[i].boundaryField()[patchI][fcI] << endl;                        
                        }
                    }
                }// end if (pfAlpha1[fcI] > 0.5)

                if(debugF_)
                {
                    Info<< "========================================================================" << nl
                        << "             Done calculating 1-phase flux for face " << faceI << nl
                        << "========================================================================" << nl << endl;
                }
            }// end if(fc_2ph_flux_needed_[faceI] && (mag(pphi[fcI]) > SMALLEST_PHI))
            faceI++;
        }// end forAll(pphi, fcI)
    }// end forAll(patches, patchI) 

    if(debugF_)
    {
        Info<< "//====================================================================\\" << nl
            << "                 Done calculating face 2-phase fluxes" << nl
            << "\\====================================================================//" << nl << nl << endl;
    }
}


void Foam::plic::calc_2ph_advFluxes
(
    const PtrList<volScalarField>& c1,
    const PtrList<volScalarField>& c0,    
    const volScalarField& h1,
    const volScalarField& h0,    
    const scalar& deltaT,
    surfaceScalarField& advFlux_rho1,
    surfaceScalarField& advFlux_rho0,
    PtrList<surfaceScalarField>& advFlux_c1,
    PtrList<surfaceScalarField>& advFlux_c0,
    surfaceScalarField& advFlux_h1,
    surfaceScalarField& advFlux_h0,
    const bool debug,
    const bool debug2,
    OFstream& os
)
{
    int n = nSpecies_;
    int i;
    label faceI, tetI, cellI, curCell_lbl;
    scalar gradAlphaOwn, gradAlphaNei, curFacePhi;
    scalar curFaceFlux_ph1, curFaceFlux_ph0, curFaceFlux_rho1, curFaceFlux_rho0, curFaceFlux_h1, curFaceFlux_h0;    
    List<scalar> curFaceFlux_c1(n-1);
    List<scalar> curFaceFlux_c0(n-1);
    scalar curTetVolFlux_ph1, curTetVolFlux_ph0, curTetFlux_rho1, curTetFlux_rho0, curTetFlux_h1, curTetFlux_h0, curTetVolSign;
    List<scalar> curTetFlux_c1(n-1);
    List<scalar> curTetFlux_c0(n-1);
    scalar curTetCurCellVolFlux_ph1, curTetCurCellVolFlux_ph0, curTetCurCellFlux_rho1, curTetCurCellFlux_rho0, curTetCurCellFlux_h1, curTetCurCellFlux_h0;
    List<scalar> curTetCurCellFlux_c1(n-1);
    List<scalar> curTetCurCellFlux_c0(n-1);

    List<tetPoints> curFaceFluxTets(20);
    List<scalar> curFaceFluxTetVols(20);
    
    const labelList& own = mesh().faceOwner();
    const labelList& nei = mesh().faceNeighbour();

    const vectorField& gradAlpha1Cells = gradAlpha1_.internalField();

    const polyBoundaryMesh& patches = mesh().boundaryMesh();

    if(debug)
    {
        os<< "//====================================================================\\" << nl
            << "                 Interpolating U field to mesh points" << nl
            << "\\====================================================================//" << nl << endl;
    }

    //U_pts_ = ptInterp_.interpolate(U_);
    ptInterp_.interpolate(U_, U_pts_);

    if(debug2)
    {
        Foam::plicFuncs::write_point_field(U_pts_, mesh());
    }

    if(debug)
    {
        os<< "//====================================================================\\" << nl
            << "                 Done interpolating U field to mesh points" << nl
            << "\\====================================================================//" << nl << endl;
    }

    if(debug)
    {
        os<< "//====================================================================\\" << nl
            << "               Collecting and distributing interface info" << nl
            << "\\====================================================================//" << nl << endl;
    }

    intfcInfo_collectData();

    c_h_collectData(c1, c0, h1, h0);

    if(debug)
    {
        os<< "//====================================================================\\" << nl
            << "            Done collecting and distributing interface info" << nl
            << "\\====================================================================//" << nl << endl;
    }

    if(debug2)
    {
        os<< "//====================================================================\\" << nl
            << "                  Tagging near-interface cells" << nl
            << "\\====================================================================//" << nl << nl
            << "------------------------------------------------------------------------" << nl
            << "cell index        mag(grad(alpha1))            near-intfc?" << nl
            << "------------------------------------------------------------------------" << endl;
    }

    forAll(cell_near_intfc_, cellI)
    {
        if(mag(gradAlpha1Cells[cellI]) > GRADALPHA_MIN)
        {
            cell_near_intfc_[cellI] = true;
        }
        else
        {
            cell_near_intfc_[cellI] = false;
        }

        if(debug2)
        {
            os<< "    " << cellI << "                  " << mag(gradAlpha1Cells[cellI]) << "                              " << cell_near_intfc_[cellI] << endl; 
        }
    }
    
    if(debug)
    {
        os<< nl
            << "//====================================================================\\" << nl
            << "                Done tagging near-interface cells" << nl
            << "\\====================================================================//" << nl << endl;
    }

    if(debug2)
    {
        os<< "//====================================================================\\" << nl
            << "                 Tagging faces requiring 2-phase flux" << nl
            << "\\====================================================================//" << nl << nl
            << "========================================================================" << nl
            << "                           Internal faces" << nl
            << "========================================================================" << nl << nl
            << "------------------------------------------------------------------------" << nl
            << "  face      own    own ni?      nei    nei ni?      need 2-phase flux?" << nl
            << "------------------------------------------------------------------------" << endl;
    }    

    for(faceI=0; faceI<mesh().nInternalFaces(); faceI++)
    {
        if(cell_near_intfc_[own[faceI]] || cell_near_intfc_[nei[faceI]])
        {
            fc_2ph_flux_needed_[faceI] = true;
        }
        else
        {
            fc_2ph_flux_needed_[faceI] = false;
        }
        
        if(debug2)
        {
            os<< "  " << faceI << "         " << own[faceI] << "      " << cell_near_intfc_[own[faceI]] << "          " << nei[faceI] << "      " << cell_near_intfc_[nei[faceI]]
                << "                  " << fc_2ph_flux_needed_[faceI] << endl;
        }
    }

    if(debug2)
    {
        os<< nl
            << "========================================================================" << nl
            << "                           Boundary faces" << nl
            << "========================================================================" << nl << endl;            
    }

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        const fvPatchVectorField& pGradAlpha1 = gradAlpha1_.boundaryField()[patchI];
        faceI = pGradAlpha1.patch().start();

        if(debug2)
        {
            os<< "--------------------------------------------------------------------" << nl
                << "                         patch: " << patchI << nl
                << "--------------------------------------------------------------------" << nl << endl;
        }

        if(pp.coupled())
        {
            if(debug2)
            {
                os<< "Coupled patch..." << nl << nl
                    << "------------------------------------------------------------------------" << nl
                    << "  face          ownGrad          neiGrad          need 2-phase flux?" << nl
                    << "------------------------------------------------------------------------" << endl;
            }

            const vectorField& gradAlphaOwn_fld = pGradAlpha1.patchInternalField();
            const vectorField& gradAlphaNei_fld = pGradAlpha1.patchNeighbourField();
            forAll(pGradAlpha1, fcI)
            {
                gradAlphaOwn = mag(gradAlphaOwn_fld[fcI]);
                gradAlphaNei = mag(gradAlphaNei_fld[fcI]);

                if(gradAlphaOwn > GRADALPHA_MIN || gradAlphaNei > GRADALPHA_MIN)
                {
                    fc_2ph_flux_needed_[faceI] = true;
                }
                else
                {
                    fc_2ph_flux_needed_[faceI] = false;
                }

                if(debug2)
                {
                    os<< "  " << faceI << "        " << gradAlphaOwn << "        " << gradAlphaNei << "        " << fc_2ph_flux_needed_[faceI] << endl;
                }

                faceI++;
            }
        }
        else if(isA<emptyPolyPatch>(pp))
        {
            if(debug2)
            {
                os<< "Empty patch..." << nl << endl;
            }

            forAll(pGradAlpha1, fcI)
            {
                fc_2ph_flux_needed_[faceI] = false;
                faceI++;
            }
        }
        else
        {
            if(debug2)
            {
                os<< "Non-coupled patch..." << nl << nl
                    << "------------------------------------------------------------------------" << nl
                    << "  face        own      own ni?          need 2-phase flux?" << nl
                    << "------------------------------------------------------------------------" << endl;
            }

            forAll(pGradAlpha1, fcI)
            {                
                if(cell_near_intfc_[own[faceI]])
                {
                    fc_2ph_flux_needed_[faceI] = true;
                }
                else
                {
                    fc_2ph_flux_needed_[faceI] = false;
                }

                if(debug2)
                {
                    os<< "  " << faceI << "           " << own[faceI] << "        " << cell_near_intfc_[own[faceI]] << "                " << fc_2ph_flux_needed_[faceI] << endl;
                }

                faceI++;
            }
        }
    }

    if(debug)
    {
        os<< "//====================================================================\\" << nl
            << "                 Done tagging faces requiring 2-phase flux" << nl
            << "\\====================================================================//" << nl << nl << endl;
    }

    if(debug)
    {
        os<< "//====================================================================\\" << nl
            << "                  Calculating face 2-phase fluxes" << nl
            << "\\====================================================================//" << nl << nl << endl;
    }    

    fAlpha1_ = fvc::interpolate(alpha_ph1_, "alpha1");
    frho1_ = fvc::interpolate(rho1_, "rho");
    frho0_ = fvc::interpolate(rho0_, "rho");
    for(i=0; i<(n-1); i++)
    {
        fc1_[i] = fvc::interpolate(c1[i], "ci");
        fc0_[i] = fvc::interpolate(c0[i], "ci");        
    }
    fh1_ = fvc::interpolate(h1, "h");
    fh0_ = fvc::interpolate(h0, "h");

    if(debug)
    {
        os<< "========================================================================" << nl
            << "                    Internal face 2-phase fluxes" << nl
            << "========================================================================" << nl << endl;
    }

    for(faceI=0; faceI<mesh().nInternalFaces(); faceI++)
    {
        curFacePhi = phi_[faceI];

        if(fc_2ph_flux_needed_[faceI] && (mag(curFacePhi) > SMALLEST_PHI))
        {
            if(debug)
            {
                os<< "========================================================================" << nl
                    << "               Calculating 2-phase flux for face " << faceI << nl
                    << "========================================================================" << nl << endl;
            }

            curFaceFlux_ph1 = 0;
            curFaceFlux_ph0 = 0;
            curFaceFlux_rho1 = 0;
            curFaceFlux_rho0 = 0;
            curFaceFlux_h1 = 0;
            curFaceFlux_h0 = 0;
            for(i=0; i<(n-1); i++)
            {
                curFaceFlux_c1[i] = 0;
                curFaceFlux_c0[i] = 0;
            }            
            
            if(debug)
            {
                os<< "========================================================================" << nl
                    << "              Step 1: Calculating face flux tetrahedrons " << nl
                    << "========================================================================" << nl << endl;
            }
            
            calcFaceFluxTets(faceI, curFacePhi, curFaceFluxTets, curFaceFluxTetVols, deltaT);

            if(debug)
            {
                os<< "========================================================================" << nl
                    << "            Done Step 1: Calculating face flux tetrahedrons" << nl
                    << "========================================================================" << nl << endl;
            }
            
            // clip each of 20 tets with the cells in the stencil for the
            // current face to find from which cell that portion of flux
            // comes from                
            const labelList& curFaceCells = faceStencil().stencil()[faceI];

            if(debug)
            {
                os<< "========================================================================" << nl
                    << "   Step 2: Intersecting face flux tetrahedrons with phase-1 subcells" << nl
                    << "========================================================================" << nl << endl;
            }

            for(tetI=0; tetI<curFaceFluxTets.size(); tetI++)
            {
                if(debug)
                {
                    os<< "--------------------------------------------------------------------" << nl
                        << "           Calculating contributions of tetrahedron " << tetI << nl
                        << "--------------------------------------------------------------------" << nl << endl;
                }

                curTetVolFlux_ph1 = 0;
                curTetVolFlux_ph0 = 0;
                curTetFlux_rho1 = 0;
                curTetFlux_rho0 = 0;
                curTetFlux_h1 = 0;
                curTetFlux_h0 = 0;
                for(i=0; i<(n-1); i++)
                {
                    curTetFlux_c1[i] = 0;
                    curTetFlux_c0[i] = 0;
                }                

                for(cellI=0; cellI<curFaceCells.size(); cellI++)
                {
                    if(debug)
                    {
                        os<< "--------------------------------------------------------------------" << nl
                            << "                   Intersecting with cell " << curFaceCells[cellI] << nl
                            << "--------------------------------------------------------------------" << nl << endl;
                    }

                    curCell_lbl = curFaceCells[cellI];
                    curTetCurCellVolFlux_ph1 = 0;
                    curTetCurCellVolFlux_ph0 = 0;
                    tet_cell_intersect(curFaceFluxTets[tetI], curCell_lbl, curTetCurCellVolFlux_ph1, curTetCurCellVolFlux_ph0);
                    curTetVolFlux_ph1 += curTetCurCellVolFlux_ph1;
                    curTetVolFlux_ph0 += curTetCurCellVolFlux_ph0;

                    curTetCurCellFlux_rho1 = rho1_flatFld_[curCell_lbl]*curTetCurCellVolFlux_ph1;
                    curTetCurCellFlux_rho0 = rho0_flatFld_[curCell_lbl]*curTetCurCellVolFlux_ph0;
                    curTetFlux_rho1 += curTetCurCellFlux_rho1;
                    curTetFlux_rho0 += curTetCurCellFlux_rho0;

                    curTetCurCellFlux_h1 = h1_flatFld_[curCell_lbl]*curTetCurCellVolFlux_ph1;
                    curTetCurCellFlux_h0 = h0_flatFld_[curCell_lbl]*curTetCurCellVolFlux_ph0;
                    curTetFlux_h1 += curTetCurCellFlux_h1;
                    curTetFlux_h0 += curTetCurCellFlux_h0;

                    if(debug)
                    {
                        os<< "curCell alpha1 = " << alpha_ph1_flatFld_[curCell_lbl] << nl
                            << "curTetCurCellVolFlux_ph1 = " << curTetCurCellVolFlux_ph1 << "  curTetCurCellVolFlux_ph0 = " << curTetCurCellVolFlux_ph0 << nl
                            << "curCell rho1 = " << rho1_flatFld_[curCell_lbl] << "  curCell rho0 = " << rho0_flatFld_[curCell_lbl] << endl;
                    }

                    for(i=0; i<(n-1); i++)
                    {
                        curTetCurCellFlux_c1[i] = c1_flatFld_[i][curCell_lbl]*curTetCurCellVolFlux_ph1;                        
                        curTetCurCellFlux_c0[i] = c0_flatFld_[i][curCell_lbl]*curTetCurCellVolFlux_ph0;
                        curTetFlux_c1[i] += curTetCurCellFlux_c1[i];
                        curTetFlux_c0[i] += curTetCurCellFlux_c0[i];
                        
                        if(debug)
                        {
                            os<< "  c1[" << i << "] = " << c1_flatFld_[i][curCell_lbl] << "  curTet curCell c1[" << i << "] flux = " << curTetCurCellFlux_c1[i] << nl
                                << "  c0[" << i << "] = " << c0_flatFld_[i][curCell_lbl] << "  curTet curCell c0[" << i << "] flux = " << curTetCurCellFlux_c0[i] << endl;
                        }
                    }// end for(i=0; i<(n-1); i++)

                    if(debug)
                    {
                        os<< nl
                            << "--------------------------------------------------------------------" << endl;                            
                    }
                }// end for(cellI=0; cellI<curFaceCells.size(); cellI++)

                curTetVolSign = sign(curFaceFluxTetVols[tetI]);
                curFaceFlux_ph1 += curTetVolFlux_ph1*curTetVolSign;
                curFaceFlux_ph0 += curTetVolFlux_ph0*curTetVolSign;
                curFaceFlux_rho1 += curTetFlux_rho1*curTetVolSign;
                curFaceFlux_rho0 += curTetFlux_rho0*curTetVolSign;
                curFaceFlux_h1 += curTetFlux_h1*curTetVolSign;
                curFaceFlux_h0 += curTetFlux_h0*curTetVolSign;

                if(debug)
                {
                    os<< "Tetrahedron " << tetI << nl
                        << "Phase-1 flux contribution:  " << curTetVolFlux_ph1 << nl
                        << "Phase-0 flux contribution:  " << curTetVolFlux_ph0 << nl
                        << "Phase-1 + Phase-0 flux:  " << curTetVolFlux_ph1 + curTetVolFlux_ph0 << nl
                        << "Tetrahedron vol:  " << curFaceFluxTetVols[tetI] << nl
                        << "Tetrahedron vol flux error: " << mag(curTetVolFlux_ph1 + curTetVolFlux_ph0 - mag(curFaceFluxTetVols[tetI])) << nl
                        << "Sign of flux:  " << curTetVolSign << nl
                        << "Cumulative phase-1 face flux:  " << curFaceFlux_ph1 << nl
                        << "Cumulative phase-0 face flux:  " << curFaceFlux_ph0 << endl;

                    if(mag(curTetVolFlux_ph1 + curTetVolFlux_ph0 - mag(curFaceFluxTetVols[tetI])) > 1E-15)
                    {
                        os<< "Tetrahedron vol flux error > 1E-15!!"
                            << endl;
                    }
                }

                for(i=0; i<(n-1); i++)
                {
                    curFaceFlux_c1[i] += curTetFlux_c1[i]*curTetVolSign;
                    curFaceFlux_c0[i] += curTetFlux_c0[i]*curTetVolSign;

                    if(debug)
                    {
                        os<< "c1[" << i << "] flux contribution:  " << curTetFlux_c1[i]*curTetVolSign << nl
                            << "Cumulative c1[" << i << "] face flux:  " << curFaceFlux_c1[i] << nl
                            << "c0[" << i << "] flux contribution:  " << curTetFlux_c0[i]*curTetVolSign << nl
                            << "Cumulative c0[" << i << "] face flux:  " << curFaceFlux_c0[i] << endl;                            
                    }
                }// end for(i=0; i<(n-1); i++)

                if(debug)
                {
                    os<< nl
                        << "--------------------------------------------------------------------" << endl;                            
                }
            }// end for(label tetI=0; tetI<curFaceFluxTets.size(); tetI++)

            if(debug)
            {
                os<< "========================================================================" << nl
                    << " Done Step 2: Intersecting face flux tetrahedrons with phase-1 subcells" << nl
                    << "========================================================================" << nl << endl;
            }

            phiAlpha1_[faceI] = curFaceFlux_ph1/deltaT;
            phiAlpha0_[faceI] = curFaceFlux_ph0/deltaT;
            advFlux_rho1[faceI] = curFaceFlux_rho1/deltaT;
            advFlux_rho0[faceI] = curFaceFlux_rho0/deltaT;
            advFlux_h1[faceI] = curFaceFlux_h1/deltaT;
            advFlux_h0[faceI] = curFaceFlux_h0/deltaT;

            if(debug)
            {
                os<< "Total phase-1 face flux:  " << phiAlpha1_[faceI] << nl
                    << "Total phase-0 face flux:  " << phiAlpha0_[faceI] << nl
                    << "Total phase-1 + phase-0 face flux:  " << phiAlpha1_[faceI] + phiAlpha0_[faceI] << nl
                    << "face phi:  " << curFacePhi << nl
                    << "face phi error:  " << mag(phiAlpha1_[faceI] + phiAlpha0_[faceI] - curFacePhi) << endl;

                if(mag(phiAlpha1_[faceI] + phiAlpha0_[faceI] - curFacePhi) > 1E-12)
                {
                    os<< "face phi error > 1E-12!!"
                        << endl;
                }
            }

            for(i=0; i<(n-1); i++)
            {
                advFlux_c1[i][faceI] = curFaceFlux_c1[i]/deltaT;
                advFlux_c0[i][faceI] = curFaceFlux_c0[i]/deltaT;

                if(debug)
                {
                    os<< "Total c1[" << i << "] face flux:  " << advFlux_c1[i][faceI] << nl
                        << "Total c0[" << i << "] face flux:  " << advFlux_c0[i][faceI] << endl;                        
                }
            }            

            if(debug)
            {
                os<< nl
                    << "========================================================================" << nl
                    << "              Done calculating 2-phase flux for face " << faceI << nl
                    << "========================================================================" << nl << endl;
            }
        }// end if(fc_2ph_flux_needed_[faceI] && (mag(phi_[faceI]) > SMALLEST_PHI))
        else
        {
            if(debug)
            {
                os<< "========================================================================" << nl
                    << "               Calculating 1-phase flux for face " << faceI << nl
                    << "========================================================================" << nl << endl;
            }

            phiAlpha1_[faceI] = phi_[faceI]*fAlpha1_[faceI];
            phiAlpha0_[faceI] = phi_[faceI]*(1 - fAlpha1_[faceI]);            

            if(debug)
            {
                os<< "Face alpha1:  " << fAlpha1_[faceI] << nl
                    << "Total phase-1 face flux:  " << phiAlpha1_[faceI] << nl
                    << "Total phase-0 face flux:  " << phiAlpha0_[faceI] << nl
                    << "Total phase-1 + phase-0 face flux:  " << phiAlpha1_[faceI] + phiAlpha0_[faceI] << nl
                    << "face phi:  " << phi_[faceI] << nl
                    << "face phi error:  " << mag(phiAlpha1_[faceI] + phiAlpha0_[faceI] - phi_[faceI]) << endl;

                if(mag(phiAlpha1_[faceI] + phiAlpha0_[faceI] - phi_[faceI]) > 1E-12)
                {
                    os<< "face phi error > 1E-12!!"
                        << endl;
                }
            }
            
            if(fAlpha1_[faceI] > 0.5)
            {
                advFlux_rho1[faceI] = phi_[faceI]*frho1_[faceI];
                advFlux_rho0[faceI] = 0;
                advFlux_h1[faceI] = phi_[faceI]*fh1_[faceI];
                advFlux_h0[faceI] = 0;

                for(i=0; i<(n-1); i++)
                {
                    advFlux_c1[i][faceI] = phi_[faceI]*fc1_[i][faceI];
                    advFlux_c0[i][faceI] = 0;

                    if(debug)
                    {
                        os<< "Total c1[" << i << "] face flux:  " << advFlux_c1[i][faceI] << nl
                            << "Total c0[" << i << "] face flux:  " << advFlux_c0[i][faceI] << endl;                        
                    }
                }
            }// end if (fAlpha1_[faceI] > 0.5)
            else
            {
                advFlux_rho1[faceI] = 0;
                advFlux_rho0[faceI] = phi_[faceI]*frho0_[faceI];
                advFlux_h1[faceI] = 0;
                advFlux_h0[faceI] = phi_[faceI]*fh0_[faceI];

                for(i=0; i<(n-1); i++)
                {
                    advFlux_c1[i][faceI] = 0;
                    advFlux_c0[i][faceI] = phi_[faceI]*fc0_[i][faceI];

                    if(debug)
                    {
                        os<< "Total c1[" << i << "] face flux:  " << advFlux_c1[i][faceI] << nl
                            << "Total c0[" << i << "] face flux:  " << advFlux_c0[i][faceI] << endl;                        
                    }
                }
            }// end if (fAlpha1_[faceI] > 0.5)

            if(debug)
            {
                os<< "========================================================================" << nl
                    << "             Done calculating 1-phase flux for face " << faceI << nl
                    << "========================================================================" << nl << endl;
            }
        }// end !if(fc_2ph_flux_needed_[faceI] && (mag(phi_[faceI]) > SMALLEST_PHI))
    }

    if(debug)
    {
        os<< nl
            << "========================================================================" << nl
            << "                    Boundary face 2-phase fluxes" << nl
            << "========================================================================" << nl << endl;
    }

    forAll(phi_.boundaryField(), patchI)
    {
        const polyPatch& pp = patches[patchI];
        const fvsPatchScalarField& pphi = phi_.boundaryField()[patchI];
        const fvsPatchScalarField& pfAlpha1 = fAlpha1_.boundaryField()[patchI];
        const fvsPatchScalarField& pfrho1 = frho1_.boundaryField()[patchI];
        const fvsPatchScalarField& pfrho0 = frho0_.boundaryField()[patchI];
        const fvsPatchScalarField& pfh1 = fh1_.boundaryField()[patchI];
        const fvsPatchScalarField& pfh0 = fh0_.boundaryField()[patchI];
        fvsPatchScalarField& pphiAlpha1 = phiAlpha1_.boundaryField()[patchI];
        fvsPatchScalarField& pphiAlpha0 = phiAlpha0_.boundaryField()[patchI];
        fvsPatchScalarField& padvFlux_rho1 = advFlux_rho1.boundaryField()[patchI];
        fvsPatchScalarField& padvFlux_rho0 = advFlux_rho0.boundaryField()[patchI];
        fvsPatchScalarField& padvFlux_h1 = advFlux_h1.boundaryField()[patchI];
        fvsPatchScalarField& padvFlux_h0 = advFlux_h0.boundaryField()[patchI];
        
        faceI = pp.start();

        if(debug)
        {
            os<< "--------------------------------------------------------------------" << nl
                << "                         patch: " << patchI << nl
                << "--------------------------------------------------------------------" << nl << endl;
        }

        forAll(pphi, fcI)
        {
            curFacePhi = pphi[fcI];

            if(fc_2ph_flux_needed_[faceI] && (mag(curFacePhi) > SMALLEST_PHI))
            {
                if(debug)
                {
                    os<< "========================================================================" << nl
                        << "               Calculating 2-phase flux for face " << faceI << nl
                        << "========================================================================" << nl << endl;
                }
                
                curFaceFlux_ph1 = 0;
                curFaceFlux_ph0 = 0;
                curFaceFlux_rho1 = 0;
                curFaceFlux_rho0 = 0;
                curFaceFlux_h1 = 0;
                curFaceFlux_h0 = 0;
                for(i=0; i<(n-1); i++)
                {
                    curFaceFlux_c1[i] = 0;
                    curFaceFlux_c0[i] = 0;
                }

                if(debug)
                {
                    os<< "========================================================================" << nl
                        << "              Step 1: Calculating face flux tetrahedrons " << nl
                        << "========================================================================" << nl << endl;
                }                            

                calcFaceFluxTets(faceI, curFacePhi, curFaceFluxTets, curFaceFluxTetVols, deltaT);
            
                if(debug)
                {
                    os<< "========================================================================" << nl
                        << "            Done Step 1: Calculating face flux tetrahedrons" << nl
                        << "========================================================================" << nl << endl;
                }

                // clip each of 20 tets with the cells in the stencil for the
                // current face to find from which cell that portion of flux
                // comes from    
                const labelList& curFaceCells = faceStencil().stencil()[faceI];

                if(debug)
                {
                    os<< "========================================================================" << nl
                        << "   Step 2: Intersecting face flux tetrahedrons with phase-1 subcells" << nl
                        << "========================================================================" << nl << endl;
                }

                for(tetI=0; tetI<curFaceFluxTets.size(); tetI++)
                {
                    if(debug)
                    {
                        os<< "--------------------------------------------------------------------" << nl
                            << "           Calculating contributions of tetrahedron " << tetI << nl
                            << "--------------------------------------------------------------------" << nl << endl;
                    }

                    curTetVolFlux_ph1 = 0;
                    curTetVolFlux_ph0 = 0;
                    curTetFlux_rho1 = 0;
                    curTetFlux_rho0 = 0;
                    curTetFlux_h1 = 0;
                    curTetFlux_h0 = 0;
                    for(i=0; i<(n-1); i++)
                    {
                        curTetFlux_c1[i] = 0;
                        curTetFlux_c0[i] = 0;
                    }

                    for(cellI=0; cellI<curFaceCells.size(); cellI++)
                    {
                        if(debug)
                        {
                            os<< "--------------------------------------------------------------------" << nl
                                << "                   Intersecting with cell " << curFaceCells[cellI] << nl
                                << "--------------------------------------------------------------------" << nl << endl;
                        }
                        
                        curCell_lbl = curFaceCells[cellI];
                        curTetCurCellVolFlux_ph1 = 0;
                        curTetCurCellVolFlux_ph0 = 0;
                        tet_cell_intersect(curFaceFluxTets[tetI], curCell_lbl, curTetCurCellVolFlux_ph1, curTetCurCellVolFlux_ph0);
                        curTetVolFlux_ph1 += curTetCurCellVolFlux_ph1;
                        curTetVolFlux_ph0 += curTetCurCellVolFlux_ph0;

                        curTetCurCellFlux_rho1 = rho1_flatFld_[curCell_lbl]*curTetCurCellVolFlux_ph1;
                        curTetCurCellFlux_rho0 = rho0_flatFld_[curCell_lbl]*curTetCurCellVolFlux_ph0;
                        curTetFlux_rho1 += curTetCurCellFlux_rho1;
                        curTetFlux_rho0 += curTetCurCellFlux_rho0;

                        curTetCurCellFlux_h1 = h1_flatFld_[curCell_lbl]*curTetCurCellVolFlux_ph1;
                        curTetCurCellFlux_h0 = h0_flatFld_[curCell_lbl]*curTetCurCellVolFlux_ph0;
                        curTetFlux_h1 += curTetCurCellFlux_h1;
                        curTetFlux_h0 += curTetCurCellFlux_h0;

                        if(debug)
                        {
                            os<< "curCell alpha1 = " << alpha_ph1_flatFld_[curCell_lbl] << nl
                                << "curTetCurCellVolFlux_ph1 = " << curTetCurCellVolFlux_ph1 << "  curTetCurCellVolFlux_ph0 = " << curTetCurCellVolFlux_ph0 << nl
                                << "curCell rho1 = " << rho1_flatFld_[curCell_lbl] << "  curCell rho0 = " << rho0_flatFld_[curCell_lbl] << endl;
                        }

                        for(i=0; i<(n-1); i++)
                        {
                            curTetCurCellFlux_c1[i] = c1_flatFld_[i][curCell_lbl]*curTetCurCellVolFlux_ph1;                        
                            curTetCurCellFlux_c0[i] = c0_flatFld_[i][curCell_lbl]*curTetCurCellVolFlux_ph0;
                            curTetFlux_c1[i] += curTetCurCellFlux_c1[i];
                            curTetFlux_c0[i] += curTetCurCellFlux_c0[i];
                        
                            if(debug)
                            {
                                os<< "  c1[" << i << "] = " << c1_flatFld_[i][curCell_lbl] << "  curTet curCell c1[" << i << "] flux = " << curTetCurCellFlux_c1[i] << nl
                                    << "  c0[" << i << "] = " << c0_flatFld_[i][curCell_lbl] << "  curTet curCell c0[" << i << "] flux = " << curTetCurCellFlux_c0[i] << endl;
                            }
                        }// end for(i=0; i<(n-1); i++)

                        if(debug)
                        {
                            os<< nl
                                << "--------------------------------------------------------------------" << endl;                            
                        }
                    }// end for(cellI=0; cellI<curFaceCells.size(); cellI++)

                    curTetVolSign = sign(curFaceFluxTetVols[tetI]);
                    curFaceFlux_ph1 += curTetVolFlux_ph1*curTetVolSign;
                    curFaceFlux_ph0 += curTetVolFlux_ph0*curTetVolSign;
                    curFaceFlux_rho1 += curTetFlux_rho1*curTetVolSign;
                    curFaceFlux_rho0 += curTetFlux_rho0*curTetVolSign;
                    curFaceFlux_h1 += curTetFlux_h1*curTetVolSign;
                    curFaceFlux_h0 += curTetFlux_h0*curTetVolSign;

                    if(debug)
                    {
                        os<< "Tetrahedron " << tetI << nl
                            << "Phase-1 flux contribution:  " << curTetVolFlux_ph1 << nl
                            << "Phase-0 flux contribution:  " << curTetVolFlux_ph0 << nl
                            << "Phase-1 + Phase-0 flux:  " << curTetVolFlux_ph1 + curTetVolFlux_ph0 << nl
                            << "Tetrahedron vol:  " << curFaceFluxTetVols[tetI] << nl
                            << "Tetrahedron vol flux error: " << mag(curTetVolFlux_ph1 + curTetVolFlux_ph0 - mag(curFaceFluxTetVols[tetI])) << nl
                            << "Sign of flux:  " << curTetVolSign << nl
                            << "Cumulative phase-1 face flux:  " << curFaceFlux_ph1 << nl
                            << "Cumulative phase-0 face flux:  " << curFaceFlux_ph0 << endl;

                        if(mag(curTetVolFlux_ph1 + curTetVolFlux_ph0 - mag(curFaceFluxTetVols[tetI])) > 1E-15)
                        {
                            os<< "Tetrahedron vol flux error > 1E-15!!"
                                << endl;
                        }
                    }

                    for(i=0; i<(n-1); i++)
                    {
                        curFaceFlux_c1[i] += curTetFlux_c1[i]*curTetVolSign;
                        curFaceFlux_c0[i] += curTetFlux_c0[i]*curTetVolSign;

                        if(debug)
                        {
                            os<< "c1[" << i << "] flux contribution:  " << curTetFlux_c1[i]*curTetVolSign << nl
                                << "Cumulative c1[" << i << "] face flux:  " << curFaceFlux_c1[i] << nl
                                << "c0[" << i << "] flux contribution:  " << curTetFlux_c0[i]*curTetVolSign << nl
                                << "Cumulative c0[" << i << "] face flux:  " << curFaceFlux_c0[i] << endl;                            
                        }
                    }// end for(i=0; i<(n-1); i++)

                    if(debug)
                    {
                        os<< nl
                            << "--------------------------------------------------------------------" << endl;                            
                    }
                }// end for(label tetI=0; tetI<curFaceFluxTets.size(); tetI++)

                if(debug)
                {
                    os<< "========================================================================" << nl
                        << " Done Step 2: Intersecting face flux tetrahedrons with phase-1 subcells" << nl
                        << "========================================================================" << nl << endl;
                }
                
                pphiAlpha1[fcI] = curFaceFlux_ph1/deltaT;
                pphiAlpha0[fcI] = curFaceFlux_ph0/deltaT;
                padvFlux_rho1[fcI] = curFaceFlux_rho1/deltaT;
                padvFlux_rho0[fcI] = curFaceFlux_rho0/deltaT;
                padvFlux_h1[fcI] = curFaceFlux_h1/deltaT;
                padvFlux_h0[fcI] = curFaceFlux_h0/deltaT;

                if(debug)
                {
                    os<< "Total phase-1 face flux:  " << pphiAlpha1[fcI] << nl
                        << "Total phase-0 face flux:  " << pphiAlpha0[fcI] << nl
                        << "phase-1 + phase-0 face flux:  " << pphiAlpha1[fcI] + pphiAlpha0[fcI] << nl
                        << "face phi:  " << curFacePhi << nl
                        << "face phi error:  " << mag(pphiAlpha1[fcI] + pphiAlpha0[fcI] - curFacePhi) << endl;

                    if(mag(pphiAlpha1[fcI] + pphiAlpha0[fcI] - curFacePhi) > 1E-12)
                    {
                        os<< "face phi error > 1E-12!!"
                            << endl;
                    }
                }

                for(i=0; i<(n-1); i++)
                {
                    advFlux_c1[i].boundaryField()[patchI][fcI] = curFaceFlux_c1[i]/deltaT;
                    advFlux_c0[i].boundaryField()[patchI][fcI] = curFaceFlux_c0[i]/deltaT;

                    if(debug)
                    {
                        os<< "Total c1[" << i << "] face flux:  " << advFlux_c1[i].boundaryField()[patchI][fcI] << nl
                            << "Total c0[" << i << "] face flux:  " << advFlux_c0[i].boundaryField()[patchI][fcI] << endl;                        
                    }
                }

                if(debug)
                {
                    os<< "========================================================================" << nl
                        << "              Done calculating 2-phase flux for face " << faceI << nl
                        << "========================================================================" << nl << endl;
                }
            }// end if(fc_2ph_flux_needed_[faceI] && (mag(pphi[fcI]) > SMALLEST_PHI))
            else
            {             
                if(debug)
                {
                    os<< "========================================================================" << nl
                        << "               Calculating 1-phase flux for face " << faceI << nl
                        << "========================================================================" << nl << endl;
                }

                pphiAlpha1[fcI] = pphi[fcI]*pfAlpha1[fcI];
                pphiAlpha0[fcI] = pphi[fcI]*(1 - pfAlpha1[fcI]);

                if(debug)
                {
                    os<< "Face alpha1:  " << pfAlpha1[fcI] << nl
                        << "Total phase-1 face flux:  " << pphiAlpha1[fcI] << nl
                        << "Total phase-0 face flux:  " << pphiAlpha0[fcI] << nl
                        << "Total phase-1 + phase-0 face flux:  " << pphiAlpha1[fcI] + pphiAlpha0[fcI] << nl
                        << "face phi:  " << pphi[fcI] << nl
                        << "face phi error:  " << mag(pphiAlpha1[fcI] + pphiAlpha0[fcI] - pphi[fcI]) << endl;

                    if(mag(pphiAlpha1[fcI] + pphiAlpha0[fcI] - pphi[fcI]) > 1E-12)
                    {
                        os<< "face phi error > 1E-12!!"
                            << endl;
                    }
                }

                if(pfAlpha1[fcI] > 0.5)
                {
                    padvFlux_rho1[fcI] = pphi[fcI]*pfrho1[fcI];
                    padvFlux_rho0[fcI] = 0;
                    padvFlux_h1[fcI] = pphi[fcI]*pfh1[fcI];
                    padvFlux_h0[fcI] = 0;

                    for(i=0; i<(n-1); i++)
                    {
                        advFlux_c1[i].boundaryField()[patchI][fcI] = pphi[fcI]*fc1_[i].boundaryField()[patchI][fcI];
                        advFlux_c0[i].boundaryField()[patchI][fcI] = 0;

                        if(debug)
                        {
                            os<< "Total c1[" << i << "] face flux:  " << advFlux_c1[i].boundaryField()[patchI][fcI] << nl
                                << "Total c0[" << i << "] face flux:  " << advFlux_c0[i].boundaryField()[patchI][fcI] << endl;                        
                        }
                    }
                }// end if (pfAlpha1[fcI] > 0.5)
                else
                {
                    padvFlux_rho1[fcI] = 0;
                    padvFlux_rho0[fcI] = pphi[fcI]*pfrho0[fcI];
                    padvFlux_h1[fcI] = 0;
                    padvFlux_h0[fcI] = pphi[fcI]*pfh0[fcI];

                    for(i=0; i<(n-1); i++)
                    {
                        advFlux_c1[i].boundaryField()[patchI][fcI] = 0;
                        advFlux_c0[i].boundaryField()[patchI][fcI] = pphi[fcI]*fc0_[i].boundaryField()[patchI][fcI];

                        if(debug)
                        {
                            os<< "Total c1[" << i << "] face flux:  " << advFlux_c1[i].boundaryField()[patchI][fcI] << nl
                                << "Total c0[" << i << "] face flux:  " << advFlux_c0[i].boundaryField()[patchI][fcI] << endl;
                        }
                    }
                }// end if (pfAlpha1[fcI] > 0.5)

                if(debug)
                {
                    os<< "========================================================================" << nl
                        << "             Done calculating 1-phase flux for face " << faceI << nl
                        << "========================================================================" << nl << endl;
                }
            }// end if(fc_2ph_flux_needed_[faceI] && (mag(pphi[fcI]) > SMALLEST_PHI))
            faceI++;
        }// end forAll(pphi, fcI)
    }// end forAll(patches, patchI) 

    if(debug)
    {
        os<< "//====================================================================\\" << nl
            << "                 Done calculating face 2-phase fluxes" << nl
            << "\\====================================================================//" << nl << nl << endl;
    }
}


void Foam::plic::calc_rAlpha0()
{
    scalar alpha1_cellI, alpha0_cellI, rAlpha0_cellI, tAlpha0;
    
    scalarField& rAlpha0Cells = rAlpha0_.internalField();
    const scalarField& alpha1Cells = alpha_ph1_.internalField();
    const scalarField& alpha0Cells = alpha_ph0_.internalField();

    forAll(rAlpha0Cells, cellI)
    {
        alpha1_cellI = alpha1Cells[cellI];
        alpha0_cellI = alpha0Cells[cellI];

        if(alpha1_cellI < ALPHA_2PH_MAX && alpha1_cellI > ALPHA_2PH_MIN)
        {
            tAlpha0 = 1.0 - alpha1_cellI;
            if(tAlpha0 < SMALL) tAlpha0 += SMALL;
            rAlpha0_cellI = alpha0_cellI/tAlpha0;
            rAlpha0_cellI = min(1.0, rAlpha0_cellI);
        }
        else
        {
            rAlpha0_cellI = 1.0;
        }

        rAlpha0Cells[cellI] = rAlpha0_cellI;
    }

    forAll(rAlpha0_.boundaryField(), patchI)
    {
        fvPatchScalarField& prAlpha0 = rAlpha0_.boundaryField()[patchI];
        const fvPatchScalarField& pAlpha1 = alpha_ph1_.boundaryField()[patchI];
        const fvPatchScalarField& pAlpha0 = alpha_ph0_.boundaryField()[patchI];

        forAll(prAlpha0, fcI)
        {
            alpha1_cellI = pAlpha1[fcI];
            alpha0_cellI = pAlpha0[fcI];

            if(alpha1_cellI < ALPHA_2PH_MAX && alpha1_cellI > ALPHA_2PH_MIN)
            {
                tAlpha0 = 1.0 - alpha1_cellI;
                if(tAlpha0 < SMALL) tAlpha0 += SMALL;
                rAlpha0_cellI = alpha0_cellI/tAlpha0;
                rAlpha0_cellI = min(1.0, rAlpha0_cellI);
            }
            else
            {
                rAlpha0_cellI = 1.0;
            }

            prAlpha0[fcI] = rAlpha0_cellI;
        }
    }
}


void Foam::plic::rAlpha0_collectData()
{
    const mapDistribute& map = faceStencil().map();

    // Insert my internal values
    forAll(rAlpha0_, cellI)
    {        
        rAlpha0_flatFld_[cellI] = rAlpha0_[cellI];
    }
    // Insert my boundary values
    forAll(rAlpha0_.boundaryField(), patchI)
    {        
        const fvPatchScalarField& prAlpha0 = rAlpha0_.boundaryField()[patchI];
        label nCompact =
            prAlpha0.patch().start()
            -rAlpha0_.mesh().nInternalFaces()
            +rAlpha0_.mesh().nCells();

        forAll(prAlpha0, faceI)
        {           
            rAlpha0_flatFld_[nCompact] = prAlpha0[faceI];
            nCompact++;
        }
    }

    // Do all swapping    
    map.distribute(rAlpha0_flatFld_);
    
    if(debug2_)
    {        
        Foam::plicFuncs::write_flatFld(rAlpha0_flatFld_, rAlpha0_);
    }
}


void Foam::plic::calc_2ph_advFluxes
(
    const PtrList<volScalarField>& c1,
    const PtrList<volScalarField>& c0,    
    const scalar& deltaT,    
    PtrList<surfaceScalarField>& advFlux_c1,
    PtrList<surfaceScalarField>& advFlux_c0,
    const bool debug,
    const bool debug2,
    OFstream& os
)
{
    int n = nSpecies_;
    int i;
    label faceI, tetI, cellI, curCell_lbl;
    scalar gradAlphaOwn, gradAlphaNei, curFacePhi;
    scalar curFaceFlux_ph1, curFaceFlux_ph0;

    List<scalar> curFaceFlux_c1(n-1);
    List<scalar> curFaceFlux_c0(n-1);
    scalar curTetVolFlux_ph1, curTetVolFlux_ph0, curTetVolSign;
    List<scalar> curTetFlux_c1(n-1);
    List<scalar> curTetFlux_c0(n-1);
    scalar curTetCurCellVolFlux_ph1, curTetCurCellVolFlux_ph0;
    List<scalar> curTetCurCellFlux_c1(n-1);
    List<scalar> curTetCurCellFlux_c0(n-1);

    List<tetPoints> curFaceFluxTets(20);
    List<scalar> curFaceFluxTetVols(20);
    
    const labelList& own = mesh().faceOwner();
    const labelList& nei = mesh().faceNeighbour();

    const vectorField& gradAlpha1Cells = gradAlpha1_.internalField();

    const polyBoundaryMesh& patches = mesh().boundaryMesh();

    if(debug)
    {
        os<< "//====================================================================\\" << nl
            << "                 Interpolating U field to mesh points" << nl
            << "\\====================================================================//" << nl << endl;
    }

    ptInterp_.interpolate(U_, U_pts_);

    if(debug2)
    {
        Foam::plicFuncs::write_point_field(U_pts_, mesh());
    }

    if(debug)
    {
        os<< "//====================================================================\\" << nl
            << "                 Done interpolating U field to mesh points" << nl
            << "\\====================================================================//" << nl << endl;
    }

    if(debug)
    {
        os<< "//====================================================================\\" << nl
            << "               Collecting and distributing interface info" << nl
            << "\\====================================================================//" << nl << endl;
    }

    calc_rAlpha0();
    rAlpha0_collectData();

    intfcInfo_collectData();

    c_collectData(c1, c0);

    if(debug)
    {
        os<< "//====================================================================\\" << nl
            << "            Done collecting and distributing interface info" << nl
            << "\\====================================================================//" << nl << endl;
    }

    if(debug2)
    {
        os<< "//====================================================================\\" << nl
            << "                  Tagging near-interface cells" << nl
            << "\\====================================================================//" << nl << nl
            << "------------------------------------------------------------------------" << nl
            << "cell index        mag(grad(alpha1))            near-intfc?" << nl
            << "------------------------------------------------------------------------" << endl;
    }

    forAll(cell_near_intfc_, cellI)
    {
        if(mag(gradAlpha1Cells[cellI]) > GRADALPHA_MIN)
        {
            cell_near_intfc_[cellI] = true;
        }
        else
        {
            cell_near_intfc_[cellI] = false;
        }

        if(debug2)
        {
            os<< "    " << cellI << "                  " << mag(gradAlpha1Cells[cellI]) << "                              " << cell_near_intfc_[cellI] << endl; 
        }
    }
    
    if(debug)
    {
        os<< nl
            << "//====================================================================\\" << nl
            << "                Done tagging near-interface cells" << nl
            << "\\====================================================================//" << nl << endl;
    }

    if(debug2)
    {
        os<< "//====================================================================\\" << nl
            << "                 Tagging faces requiring 2-phase flux" << nl
            << "\\====================================================================//" << nl << nl
            << "========================================================================" << nl
            << "                           Internal faces" << nl
            << "========================================================================" << nl << nl
            << "------------------------------------------------------------------------" << nl
            << "  face      own    own ni?      nei    nei ni?      need 2-phase flux?" << nl
            << "------------------------------------------------------------------------" << endl;
    }    

    for(faceI=0; faceI<mesh().nInternalFaces(); faceI++)
    {
        if(cell_near_intfc_[own[faceI]] || cell_near_intfc_[nei[faceI]])
        {
            fc_2ph_flux_needed_[faceI] = true;
        }
        else
        {
            fc_2ph_flux_needed_[faceI] = false;
        }
        
        if(debug2)
        {
            os<< "  " << faceI << "         " << own[faceI] << "      " << cell_near_intfc_[own[faceI]] << "          " << nei[faceI] << "      " << cell_near_intfc_[nei[faceI]]
                << "                  " << fc_2ph_flux_needed_[faceI] << endl;
        }
    }

    if(debug2)
    {
        os<< nl
            << "========================================================================" << nl
            << "                           Boundary faces" << nl
            << "========================================================================" << nl << endl;            
    }

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        const fvPatchVectorField& pGradAlpha1 = gradAlpha1_.boundaryField()[patchI];
        faceI = pGradAlpha1.patch().start();

        if(debug2)
        {
            os<< "--------------------------------------------------------------------" << nl
                << "                         patch: " << patchI << nl
                << "--------------------------------------------------------------------" << nl << endl;
        }

        if(pp.coupled())
        {
            if(debug2)
            {
                os<< "Coupled patch..." << nl << nl
                    << "------------------------------------------------------------------------" << nl
                    << "  face          ownGrad          neiGrad          need 2-phase flux?" << nl
                    << "------------------------------------------------------------------------" << endl;
            }

            const vectorField& gradAlphaOwn_fld = pGradAlpha1.patchInternalField();
            const vectorField& gradAlphaNei_fld = pGradAlpha1.patchNeighbourField();
            forAll(pGradAlpha1, fcI)
            {
                gradAlphaOwn = mag(gradAlphaOwn_fld[fcI]);
                gradAlphaNei = mag(gradAlphaNei_fld[fcI]);

                if(gradAlphaOwn > GRADALPHA_MIN || gradAlphaNei > GRADALPHA_MIN)
                {
                    fc_2ph_flux_needed_[faceI] = true;
                }
                else
                {
                    fc_2ph_flux_needed_[faceI] = false;
                }

                if(debug2)
                {
                    os<< "  " << faceI << "        " << gradAlphaOwn << "        " << gradAlphaNei << "        " << fc_2ph_flux_needed_[faceI] << endl;
                }

                faceI++;
            }
        }
        else if(isA<emptyPolyPatch>(pp))
        {
            if(debug2)
            {
                os<< "Empty patch..." << nl << endl;
            }

            forAll(pGradAlpha1, fcI)
            {
                fc_2ph_flux_needed_[faceI] = false;
                faceI++;
            }
        }
        else
        {
            if(debug2)
            {
                os<< "Non-coupled patch..." << nl << nl
                    << "------------------------------------------------------------------------" << nl
                    << "  face        own      own ni?          need 2-phase flux?" << nl
                    << "------------------------------------------------------------------------" << endl;
            }

            forAll(pGradAlpha1, fcI)
            {                
                if(cell_near_intfc_[own[faceI]])
                {
                    fc_2ph_flux_needed_[faceI] = true;
                }
                else
                {
                    fc_2ph_flux_needed_[faceI] = false;
                }

                if(debug2)
                {
                    os<< "  " << faceI << "           " << own[faceI] << "        " << cell_near_intfc_[own[faceI]] << "                " << fc_2ph_flux_needed_[faceI] << endl;
                }

                faceI++;
            }
        }
    }

    if(debug)
    {
        os<< "//====================================================================\\" << nl
            << "                 Done tagging faces requiring 2-phase flux" << nl
            << "\\====================================================================//" << nl << nl << endl;
    }

    if(debug)
    {
        os<< "//====================================================================\\" << nl
            << "                  Calculating face 2-phase fluxes" << nl
            << "\\====================================================================//" << nl << nl << endl;
    }    

    fAlpha1_ = fvc::interpolate(alpha_ph1_, "alpha1");    
    for(i=0; i<(n-1); i++)
    {
        fc1_[i] = fvc::interpolate(c1[i], "ci");
        fc0_[i] = fvc::interpolate(c0[i], "ci");        
    }

    if(debug)
    {
        os<< "========================================================================" << nl
            << "                    Internal face 2-phase fluxes" << nl
            << "========================================================================" << nl << endl;
    }

    for(faceI=0; faceI<mesh().nInternalFaces(); faceI++)
    {
        curFacePhi = phi_[faceI];

        if(fc_2ph_flux_needed_[faceI] && (mag(curFacePhi) > SMALLEST_PHI))
        {
            if(debug)
            {
                os<< "========================================================================" << nl
                    << "               Calculating 2-phase flux for face " << faceI << nl
                    << "========================================================================" << nl << endl;
            }

            curFaceFlux_ph1 = 0;
            curFaceFlux_ph0 = 0;            
            for(i=0; i<(n-1); i++)
            {
                curFaceFlux_c1[i] = 0;
                curFaceFlux_c0[i] = 0;
            }            
            
            if(debug)
            {
                os<< "========================================================================" << nl
                    << "              Step 1: Calculating face flux tetrahedrons " << nl
                    << "========================================================================" << nl << endl;
            }
            
            calcFaceFluxTets(faceI, curFacePhi, curFaceFluxTets, curFaceFluxTetVols, deltaT);

            if(debug)
            {
                os<< "========================================================================" << nl
                    << "            Done Step 1: Calculating face flux tetrahedrons" << nl
                    << "========================================================================" << nl << endl;
            }
            
            // clip each of 20 tets with the cells in the stencil for the
            // current face to find from which cell that portion of flux
            // comes from                
            const labelList& curFaceCells = faceStencil().stencil()[faceI];

            if(debug)
            {
                os<< "========================================================================" << nl
                    << "   Step 2: Intersecting face flux tetrahedrons with phase-1 subcells" << nl
                    << "========================================================================" << nl << endl;
            }

            for(tetI=0; tetI<curFaceFluxTets.size(); tetI++)
            {
                if(debug)
                {
                    os<< "--------------------------------------------------------------------" << nl
                        << "           Calculating contributions of tetrahedron " << tetI << nl
                        << "--------------------------------------------------------------------" << nl << endl;
                }

                curTetVolFlux_ph1 = 0;
                curTetVolFlux_ph0 = 0;                
                for(i=0; i<(n-1); i++)
                {
                    curTetFlux_c1[i] = 0;
                    curTetFlux_c0[i] = 0;
                }                

                for(cellI=0; cellI<curFaceCells.size(); cellI++)
                {
                    if(debug)
                    {
                        os<< "--------------------------------------------------------------------" << nl
                            << "                   Intersecting with cell " << curFaceCells[cellI] << nl
                            << "--------------------------------------------------------------------" << nl << endl;
                    }

                    curCell_lbl = curFaceCells[cellI];
                    curTetCurCellVolFlux_ph1 = 0;
                    curTetCurCellVolFlux_ph0 = 0;
                    tet_cell_intersect(curFaceFluxTets[tetI], curCell_lbl, curTetCurCellVolFlux_ph1, curTetCurCellVolFlux_ph0);
                    curTetVolFlux_ph1 += curTetCurCellVolFlux_ph1;
                    curTetVolFlux_ph0 += curTetCurCellVolFlux_ph0;
                    
                    if(debug)
                    {
                        os<< "curCell alpha1 = " << alpha_ph1_flatFld_[curCell_lbl] << nl
                            << "curTetCurCellVolFlux_ph1 = " << curTetCurCellVolFlux_ph1 << "  curTetCurCellVolFlux_ph0 = " << curTetCurCellVolFlux_ph0 << nl
                            << "curCell rho1 = " << rho1_flatFld_[curCell_lbl] << "  curCell rho0 = " << rho0_flatFld_[curCell_lbl] << endl;
                    }

                    for(i=0; i<(n-1); i++)
                    {
                        curTetCurCellFlux_c1[i] = c1_flatFld_[i][curCell_lbl]*curTetCurCellVolFlux_ph1;                        
                        curTetCurCellFlux_c0[i] = rAlpha0_flatFld_[curCell_lbl]*c0_flatFld_[i][curCell_lbl]*curTetCurCellVolFlux_ph0;
                        curTetFlux_c1[i] += curTetCurCellFlux_c1[i];
                        curTetFlux_c0[i] += curTetCurCellFlux_c0[i];
                        
                        if(debug)
                        {
                            os<< "  c1[" << i << "] = " << c1_flatFld_[i][curCell_lbl] << "  curTet curCell c1[" << i << "] flux = " << curTetCurCellFlux_c1[i] << nl
                                << "  c0[" << i << "] = " << c0_flatFld_[i][curCell_lbl] << "  curTet curCell c0[" << i << "] flux = " << curTetCurCellFlux_c0[i] << endl;
                        }
                    }// end for(i=0; i<(n-1); i++)

                    if(debug)
                    {
                        os<< nl
                            << "--------------------------------------------------------------------" << endl;                            
                    }
                }// end for(cellI=0; cellI<curFaceCells.size(); cellI++)

                curTetVolSign = sign(curFaceFluxTetVols[tetI]);
                curFaceFlux_ph1 += curTetVolFlux_ph1*curTetVolSign;
                curFaceFlux_ph0 += curTetVolFlux_ph0*curTetVolSign;                

                if(debug)
                {
                    os<< "Tetrahedron " << tetI << nl
                        << "Phase-1 flux contribution:  " << curTetVolFlux_ph1 << nl
                        << "Phase-0 flux contribution:  " << curTetVolFlux_ph0 << nl
                        << "Phase-1 + Phase-0 flux:  " << curTetVolFlux_ph1 + curTetVolFlux_ph0 << nl
                        << "Tetrahedron vol:  " << curFaceFluxTetVols[tetI] << nl
                        << "Tetrahedron vol flux error: " << mag(curTetVolFlux_ph1 + curTetVolFlux_ph0 - mag(curFaceFluxTetVols[tetI])) << nl
                        << "Sign of flux:  " << curTetVolSign << nl
                        << "Cumulative phase-1 face flux:  " << curFaceFlux_ph1 << nl
                        << "Cumulative phase-0 face flux:  " << curFaceFlux_ph0 << endl;

                    if(mag(curTetVolFlux_ph1 + curTetVolFlux_ph0 - mag(curFaceFluxTetVols[tetI])) > 1E-15)
                    {
                        os<< "Tetrahedron vol flux error > 1E-15!!"
                            << endl;
                    }
                }

                for(i=0; i<(n-1); i++)
                {
                    curFaceFlux_c1[i] += curTetFlux_c1[i]*curTetVolSign;
                    curFaceFlux_c0[i] += curTetFlux_c0[i]*curTetVolSign;

                    if(debug)
                    {
                        os<< "c1[" << i << "] flux contribution:  " << curTetFlux_c1[i]*curTetVolSign << nl
                            << "Cumulative c1[" << i << "] face flux:  " << curFaceFlux_c1[i] << nl
                            << "c0[" << i << "] flux contribution:  " << curTetFlux_c0[i]*curTetVolSign << nl
                            << "Cumulative c0[" << i << "] face flux:  " << curFaceFlux_c0[i] << endl;                            
                    }
                }// end for(i=0; i<(n-1); i++)

                if(debug)
                {
                    os<< nl
                        << "--------------------------------------------------------------------" << endl;                            
                }
            }// end for(label tetI=0; tetI<curFaceFluxTets.size(); tetI++)

            if(debug)
            {
                os<< "========================================================================" << nl
                    << " Done Step 2: Intersecting face flux tetrahedrons with phase-1 subcells" << nl
                    << "========================================================================" << nl << endl;
            }

            phiAlpha1_[faceI] = curFaceFlux_ph1/deltaT;
            phiAlpha0_[faceI] = curFaceFlux_ph0/deltaT;

            if(debug)
            {
                os<< "Total phase-1 face flux:  " << phiAlpha1_[faceI] << nl
                    << "Total phase-0 face flux:  " << phiAlpha0_[faceI] << nl
                    << "Total phase-1 + phase-0 face flux:  " << phiAlpha1_[faceI] + phiAlpha0_[faceI] << nl
                    << "face phi:  " << curFacePhi << nl
                    << "face phi error:  " << mag(phiAlpha1_[faceI] + phiAlpha0_[faceI] - curFacePhi) << endl;

                if(mag(phiAlpha1_[faceI] + phiAlpha0_[faceI] - curFacePhi) > 1E-12)
                {
                    os<< "face phi error > 1E-12!!"
                        << endl;
                }
            }

            for(i=0; i<(n-1); i++)
            {
                advFlux_c1[i][faceI] = curFaceFlux_c1[i]/deltaT;
                advFlux_c0[i][faceI] = curFaceFlux_c0[i]/deltaT;

                if(debug)
                {
                    os<< "Total c1[" << i << "] face flux:  " << advFlux_c1[i][faceI] << nl
                        << "Total c0[" << i << "] face flux:  " << advFlux_c0[i][faceI] << endl;                        
                }
            }            

            if(debug)
            {
                os<< nl
                    << "========================================================================" << nl
                    << "              Done calculating 2-phase flux for face " << faceI << nl
                    << "========================================================================" << nl << endl;
            }
        }// end if(fc_2ph_flux_needed_[faceI] && (mag(phi_[faceI]) > SMALLEST_PHI))
        else
        {
            if(debug)
            {
                os<< "========================================================================" << nl
                    << "               Calculating 1-phase flux for face " << faceI << nl
                    << "========================================================================" << nl << endl;
            }

            phiAlpha1_[faceI] = phi_[faceI]*fAlpha1_[faceI];
            phiAlpha0_[faceI] = phi_[faceI]*(1 - fAlpha1_[faceI]);            

            if(debug)
            {
                os<< "Face alpha1:  " << fAlpha1_[faceI] << nl
                    << "Total phase-1 face flux:  " << phiAlpha1_[faceI] << nl
                    << "Total phase-0 face flux:  " << phiAlpha0_[faceI] << nl
                    << "Total phase-1 + phase-0 face flux:  " << phiAlpha1_[faceI] + phiAlpha0_[faceI] << nl
                    << "face phi:  " << phi_[faceI] << nl
                    << "face phi error:  " << mag(phiAlpha1_[faceI] + phiAlpha0_[faceI] - phi_[faceI]) << endl;

                if(mag(phiAlpha1_[faceI] + phiAlpha0_[faceI] - phi_[faceI]) > 1E-12)
                {
                    os<< "face phi error > 1E-12!!"
                        << endl;
                }
            }
            
            if(fAlpha1_[faceI] > 0.5)
            {
                for(i=0; i<(n-1); i++)
                {
                    advFlux_c1[i][faceI] = phi_[faceI]*fc1_[i][faceI];
                    advFlux_c0[i][faceI] = 0;

                    if(debug)
                    {
                        os<< "Total c1[" << i << "] face flux:  " << advFlux_c1[i][faceI] << nl
                            << "Total c0[" << i << "] face flux:  " << advFlux_c0[i][faceI] << endl;                        
                    }
                }
            }// end if (fAlpha1_[faceI] > 0.5)
            else
            {
                for(i=0; i<(n-1); i++)
                {
                    advFlux_c1[i][faceI] = 0;
                    advFlux_c0[i][faceI] = phi_[faceI]*fc0_[i][faceI];

                    if(debug)
                    {
                        os<< "Total c1[" << i << "] face flux:  " << advFlux_c1[i][faceI] << nl
                            << "Total c0[" << i << "] face flux:  " << advFlux_c0[i][faceI] << endl;                        
                    }
                }
            }// end if (fAlpha1_[faceI] > 0.5)

            if(debug)
            {
                os<< "========================================================================" << nl
                    << "             Done calculating 1-phase flux for face " << faceI << nl
                    << "========================================================================" << nl << endl;
            }
        }// end !if(fc_2ph_flux_needed_[faceI] && (mag(phi_[faceI]) > SMALLEST_PHI))
    }

    if(debug)
    {
        os<< nl
            << "========================================================================" << nl
            << "                    Boundary face 2-phase fluxes" << nl
            << "========================================================================" << nl << endl;
    }

    forAll(phi_.boundaryField(), patchI)
    {
        const polyPatch& pp = patches[patchI];
        const fvsPatchScalarField& pphi = phi_.boundaryField()[patchI];
        const fvsPatchScalarField& pfAlpha1 = fAlpha1_.boundaryField()[patchI];        
        fvsPatchScalarField& pphiAlpha1 = phiAlpha1_.boundaryField()[patchI];
        fvsPatchScalarField& pphiAlpha0 = phiAlpha0_.boundaryField()[patchI];
        
        faceI = pp.start();

        if(debug)
        {
            os<< "--------------------------------------------------------------------" << nl
                << "                         patch: " << patchI << nl
                << "--------------------------------------------------------------------" << nl << endl;
        }

        forAll(pphi, fcI)
        {
            curFacePhi = pphi[fcI];

            if(fc_2ph_flux_needed_[faceI] && (mag(curFacePhi) > SMALLEST_PHI))
            {
                if(debug)
                {
                    os<< "========================================================================" << nl
                        << "               Calculating 2-phase flux for face " << faceI << nl
                        << "========================================================================" << nl << endl;
                }
                
                curFaceFlux_ph1 = 0;
                curFaceFlux_ph0 = 0;
                for(i=0; i<(n-1); i++)
                {
                    curFaceFlux_c1[i] = 0;
                    curFaceFlux_c0[i] = 0;
                }

                if(debug)
                {
                    os<< "========================================================================" << nl
                        << "              Step 1: Calculating face flux tetrahedrons " << nl
                        << "========================================================================" << nl << endl;
                }                            

                calcFaceFluxTets(faceI, curFacePhi, curFaceFluxTets, curFaceFluxTetVols, deltaT);
            
                if(debug)
                {
                    os<< "========================================================================" << nl
                        << "            Done Step 1: Calculating face flux tetrahedrons" << nl
                        << "========================================================================" << nl << endl;
                }

                // clip each of 20 tets with the cells in the stencil for the
                // current face to find from which cell that portion of flux
                // comes from    
                const labelList& curFaceCells = faceStencil().stencil()[faceI];

                if(debug)
                {
                    os<< "========================================================================" << nl
                        << "   Step 2: Intersecting face flux tetrahedrons with phase-1 subcells" << nl
                        << "========================================================================" << nl << endl;
                }

                for(tetI=0; tetI<curFaceFluxTets.size(); tetI++)
                {
                    if(debug)
                    {
                        os<< "--------------------------------------------------------------------" << nl
                            << "           Calculating contributions of tetrahedron " << tetI << nl
                            << "--------------------------------------------------------------------" << nl << endl;
                    }

                    curTetVolFlux_ph1 = 0;
                    curTetVolFlux_ph0 = 0;
                    for(i=0; i<(n-1); i++)
                    {
                        curTetFlux_c1[i] = 0;
                        curTetFlux_c0[i] = 0;
                    }

                    for(cellI=0; cellI<curFaceCells.size(); cellI++)
                    {
                        if(debug)
                        {
                            os<< "--------------------------------------------------------------------" << nl
                                << "                   Intersecting with cell " << curFaceCells[cellI] << nl
                                << "--------------------------------------------------------------------" << nl << endl;
                        }
                        
                        curCell_lbl = curFaceCells[cellI];
                        curTetCurCellVolFlux_ph1 = 0;
                        curTetCurCellVolFlux_ph0 = 0;
                        tet_cell_intersect(curFaceFluxTets[tetI], curCell_lbl, curTetCurCellVolFlux_ph1, curTetCurCellVolFlux_ph0);
                        curTetVolFlux_ph1 += curTetCurCellVolFlux_ph1;
                        curTetVolFlux_ph0 += curTetCurCellVolFlux_ph0;                        

                        if(debug)
                        {
                            os<< "curCell alpha1 = " << alpha_ph1_flatFld_[curCell_lbl] << nl
                                << "curTetCurCellVolFlux_ph1 = " << curTetCurCellVolFlux_ph1 << "  curTetCurCellVolFlux_ph0 = " << curTetCurCellVolFlux_ph0 << nl
                                << "curCell rho1 = " << rho1_flatFld_[curCell_lbl] << "  curCell rho0 = " << rho0_flatFld_[curCell_lbl] << endl;
                        }

                        for(i=0; i<(n-1); i++)
                        {
                            curTetCurCellFlux_c1[i] = c1_flatFld_[i][curCell_lbl]*curTetCurCellVolFlux_ph1;                        
                            curTetCurCellFlux_c0[i] = rAlpha0_flatFld_[curCell_lbl]*c0_flatFld_[i][curCell_lbl]*curTetCurCellVolFlux_ph0;
                            curTetFlux_c1[i] += curTetCurCellFlux_c1[i];
                            curTetFlux_c0[i] += curTetCurCellFlux_c0[i];
                        
                            if(debug)
                            {
                                os<< "  c1[" << i << "] = " << c1_flatFld_[i][curCell_lbl] << "  curTet curCell c1[" << i << "] flux = " << curTetCurCellFlux_c1[i] << nl
                                    << "  c0[" << i << "] = " << c0_flatFld_[i][curCell_lbl] << "  curTet curCell c0[" << i << "] flux = " << curTetCurCellFlux_c0[i] << endl;
                            }
                        }// end for(i=0; i<(n-1); i++)

                        if(debug)
                        {
                            os<< nl
                                << "--------------------------------------------------------------------" << endl;                            
                        }
                    }// end for(cellI=0; cellI<curFaceCells.size(); cellI++)

                    curTetVolSign = sign(curFaceFluxTetVols[tetI]);
                    curFaceFlux_ph1 += curTetVolFlux_ph1*curTetVolSign;
                    curFaceFlux_ph0 += curTetVolFlux_ph0*curTetVolSign;

                    if(debug)
                    {
                        os<< "Tetrahedron " << tetI << nl
                            << "Phase-1 flux contribution:  " << curTetVolFlux_ph1 << nl
                            << "Phase-0 flux contribution:  " << curTetVolFlux_ph0 << nl
                            << "Phase-1 + Phase-0 flux:  " << curTetVolFlux_ph1 + curTetVolFlux_ph0 << nl
                            << "Tetrahedron vol:  " << curFaceFluxTetVols[tetI] << nl
                            << "Tetrahedron vol flux error: " << mag(curTetVolFlux_ph1 + curTetVolFlux_ph0 - mag(curFaceFluxTetVols[tetI])) << nl
                            << "Sign of flux:  " << curTetVolSign << nl
                            << "Cumulative phase-1 face flux:  " << curFaceFlux_ph1 << nl
                            << "Cumulative phase-0 face flux:  " << curFaceFlux_ph0 << endl;

                        if(mag(curTetVolFlux_ph1 + curTetVolFlux_ph0 - mag(curFaceFluxTetVols[tetI])) > 1E-15)
                        {
                            os<< "Tetrahedron vol flux error > 1E-15!!"
                                << endl;
                        }
                    }

                    for(i=0; i<(n-1); i++)
                    {
                        curFaceFlux_c1[i] += curTetFlux_c1[i]*curTetVolSign;
                        curFaceFlux_c0[i] += curTetFlux_c0[i]*curTetVolSign;

                        if(debug)
                        {
                            os<< "c1[" << i << "] flux contribution:  " << curTetFlux_c1[i]*curTetVolSign << nl
                                << "Cumulative c1[" << i << "] face flux:  " << curFaceFlux_c1[i] << nl
                                << "c0[" << i << "] flux contribution:  " << curTetFlux_c0[i]*curTetVolSign << nl
                                << "Cumulative c0[" << i << "] face flux:  " << curFaceFlux_c0[i] << endl;                            
                        }
                    }// end for(i=0; i<(n-1); i++)

                    if(debug)
                    {
                        os<< nl
                            << "--------------------------------------------------------------------" << endl;                            
                    }
                }// end for(label tetI=0; tetI<curFaceFluxTets.size(); tetI++)

                if(debug)
                {
                    os<< "========================================================================" << nl
                        << " Done Step 2: Intersecting face flux tetrahedrons with phase-1 subcells" << nl
                        << "========================================================================" << nl << endl;
                }
                
                pphiAlpha1[fcI] = curFaceFlux_ph1/deltaT;
                pphiAlpha0[fcI] = curFaceFlux_ph0/deltaT;

                if(debug)
                {
                    os<< "Total phase-1 face flux:  " << pphiAlpha1[fcI] << nl
                        << "Total phase-0 face flux:  " << pphiAlpha0[fcI] << nl
                        << "phase-1 + phase-0 face flux:  " << pphiAlpha1[fcI] + pphiAlpha0[fcI] << nl
                        << "face phi:  " << curFacePhi << nl
                        << "face phi error:  " << mag(pphiAlpha1[fcI] + pphiAlpha0[fcI] - curFacePhi) << endl;

                    if(mag(pphiAlpha1[fcI] + pphiAlpha0[fcI] - curFacePhi) > 1E-12)
                    {
                        os<< "face phi error > 1E-12!!"
                            << endl;
                    }
                }

                for(i=0; i<(n-1); i++)
                {
                    advFlux_c1[i].boundaryField()[patchI][fcI] = curFaceFlux_c1[i]/deltaT;
                    advFlux_c0[i].boundaryField()[patchI][fcI] = curFaceFlux_c0[i]/deltaT;

                    if(debug)
                    {
                        os<< "Total c1[" << i << "] face flux:  " << advFlux_c1[i].boundaryField()[patchI][fcI] << nl
                            << "Total c0[" << i << "] face flux:  " << advFlux_c0[i].boundaryField()[patchI][fcI] << endl;                        
                    }
                }

                if(debug)
                {
                    os<< "========================================================================" << nl
                        << "              Done calculating 2-phase flux for face " << faceI << nl
                        << "========================================================================" << nl << endl;
                }
            }// end if(fc_2ph_flux_needed_[faceI] && (mag(pphi[fcI]) > SMALLEST_PHI))
            else
            {             
                if(debug)
                {
                    os<< "========================================================================" << nl
                        << "               Calculating 1-phase flux for face " << faceI << nl
                        << "========================================================================" << nl << endl;
                }

                pphiAlpha1[fcI] = pphi[fcI]*pfAlpha1[fcI];
                pphiAlpha0[fcI] = pphi[fcI]*(1 - pfAlpha1[fcI]);

                if(debug)
                {
                    os<< "Face alpha1:  " << pfAlpha1[fcI] << nl
                        << "Total phase-1 face flux:  " << pphiAlpha1[fcI] << nl
                        << "Total phase-0 face flux:  " << pphiAlpha0[fcI] << nl
                        << "Total phase-1 + phase-0 face flux:  " << pphiAlpha1[fcI] + pphiAlpha0[fcI] << nl
                        << "face phi:  " << pphi[fcI] << nl
                        << "face phi error:  " << mag(pphiAlpha1[fcI] + pphiAlpha0[fcI] - pphi[fcI]) << endl;

                    if(mag(pphiAlpha1[fcI] + pphiAlpha0[fcI] - pphi[fcI]) > 1E-12)
                    {
                        os<< "face phi error > 1E-12!!"
                            << endl;
                    }
                }

                if(pfAlpha1[fcI] > 0.5)
                {                    
                    for(i=0; i<(n-1); i++)
                    {
                        advFlux_c1[i].boundaryField()[patchI][fcI] = pphi[fcI]*fc1_[i].boundaryField()[patchI][fcI];
                        advFlux_c0[i].boundaryField()[patchI][fcI] = 0;

                        if(debug)
                        {
                            os<< "Total c1[" << i << "] face flux:  " << advFlux_c1[i].boundaryField()[patchI][fcI] << nl
                                << "Total c0[" << i << "] face flux:  " << advFlux_c0[i].boundaryField()[patchI][fcI] << endl;                        
                        }
                    }
                }// end if (pfAlpha1[fcI] > 0.5)
                else
                {
                    for(label i=0; i<(nSpecies_ - 1); i++)
                    {
                        advFlux_c1[i].boundaryField()[patchI][fcI] = 0;
                        advFlux_c0[i].boundaryField()[patchI][fcI] = pphi[fcI]*fc0_[i].boundaryField()[patchI][fcI];

                        if(debug)
                        {
                            os<< "Total c1[" << i << "] face flux:  " << advFlux_c1[i].boundaryField()[patchI][fcI] << nl
                                << "Total c0[" << i << "] face flux:  " << advFlux_c0[i].boundaryField()[patchI][fcI] << endl;
                        }
                    }
                }// end if (pfAlpha1[fcI] > 0.5)

                if(debug)
                {
                    os<< "========================================================================" << nl
                        << "             Done calculating 1-phase flux for face " << faceI << nl
                        << "========================================================================" << nl << endl;
                }
            }// end if(fc_2ph_flux_needed_[faceI] && (mag(pphi[fcI]) > SMALLEST_PHI))
            faceI++;
        }// end forAll(pphi, fcI)
    }// end forAll(patches, patchI) 

    if(debug)
    {
        os<< "//====================================================================\\" << nl
            << "                 Done calculating face 2-phase fluxes" << nl
            << "\\====================================================================//" << nl << nl << endl;
    }
}


void Foam::plic::plic_correct()
{
    intfc_correct();
    calc_face_phaseFluxes();
}


bool Foam::plic::isEq_pt(const point& pt1, const point& pt2)
{
    if((mag(pt1.x()-pt2.x()) < PT_DIST_TOL) && (mag(pt1.y()-pt2.y()) < PT_DIST_TOL) && (mag(pt1.z()-pt2.z()) < PT_DIST_TOL))
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
