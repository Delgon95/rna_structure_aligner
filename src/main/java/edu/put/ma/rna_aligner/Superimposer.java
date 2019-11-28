package edu.put.ma.rna_aligner;

import java.util.ArrayList;

import org.biojava.nbio.structure.jama.Matrix;
import org.biojava.nbio.structure.jama.SingularValueDecomposition;

// Rewritten SVDSuperimposer class from BioJava using just 'Coordinates' instead of 'Atom' structure
// Also inputs ArrayList instead of Array.
public class Superimposer {
    Matrix rot;
    Matrix tran;
    
    Coordinates translation;
    Coordinates centroidA;
    Coordinates centroidB;

    public Superimposer(final ArrayList<Coordinates> atomSet1, final ArrayList<Coordinates> atomSet2) {
    	centroidA = Coordinates.getCentroid(atomSet1);
    	centroidB = Coordinates.getCentroid(atomSet2);
    	
    	// center at centroid
    	Matrix centereda = CenteredCoordsMatrix(atomSet1, centroidA);
    	Matrix centeredb = CenteredCoordsMatrixTransposed(atomSet2, centroidB);    	
    	
    	calculate(centereda, centeredb);
    }
    
    private final Matrix CenteredCoordsMatrix(final ArrayList<Coordinates> atomSet, final Coordinates centroid) {
    	Coordinates shiftVector = Coordinates.getCenterVector(atomSet, centroid);
    	Matrix centered = new Matrix(atomSet.size(), 3);
    	
    	for (int i = 0; i < atomSet.size(); ++i) {
    		centered.set(i, 0, atomSet.get(i).x + shiftVector.x);
    		centered.set(i, 1, atomSet.get(i).y + shiftVector.y);
    		centered.set(i, 2, atomSet.get(i).z + shiftVector.z);
    	}
    	return centered;
    }
    
    private final Matrix CenteredCoordsMatrixTransposed(final ArrayList<Coordinates> atomSet, final Coordinates centroid) {
    	Coordinates shiftVector = Coordinates.getCenterVector(atomSet, centroid);
    	Matrix centered = new Matrix(3, atomSet.size());
    	
    	for (int i = 0; i < atomSet.size(); ++i) {
    		centered.set(0, i, atomSet.get(i).x + shiftVector.x);
    		centered.set(1, i, atomSet.get(i).y + shiftVector.y);
    		centered.set(2, i, atomSet.get(i).z + shiftVector.z);
    	}
    	return centered;
    }

    private void calculate (Matrix a, Matrix b_trans) {
    	Matrix corr = b_trans.times(a);

    	SingularValueDecomposition svd = corr.svd();

    	Matrix u = svd.getU();
    	// v is alreaady transposed ! difference to numermic python ...
    	Matrix vt = svd.getV();

    	Matrix u_transp = u.transpose();

    	Matrix rot_nottrans = vt.times(u_transp);
    	rot = rot_nottrans.transpose();

    	// check if we have found a reflection

    	double det = rot.det();

    	if (det < 0) {
    		vt = vt.transpose();
    		vt.set(2, 0, (0 - vt.get(2, 0)));
    		vt.set(2, 1, (0 - vt.get(2, 1)));
    		vt.set(2, 2, (0 - vt.get(2, 2)));

    		Matrix nv_transp = vt.transpose();
    		rot_nottrans = nv_transp.times(u_transp);
    		rot = rot_nottrans.transpose();

    	}

    	Coordinates cb_tmp = centroidB.timesMatrix(rot);
    	translation = centroidA.clone().substract(cb_tmp);
    }
    
    public final Matrix getRotation() {
    	return rot;
    }
    
    public final Coordinates getTranslation(){
    	return translation;
    }
    
    public static double getRMS(final ArrayList<Coordinates> atomSet1, final ArrayList<Coordinates> atomSet2) {
    	double sum = 0.0;
    	for (int i = 0 ; i < atomSet1.size(); i++) {
    		double d = Coordinates.getDistance(atomSet1.get(i), atomSet2.get(i));
    		sum += (d*d);
    	}

    	double avd = (sum / atomSet1.size());
    	return Math.sqrt(avd);
    }
}
