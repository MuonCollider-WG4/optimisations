//	BLCMDsolenoid.cc
/*
This source file is part of G4beamline, http://g4beamline.muonsinc.com
Copyright (C) 2002-2013 by Tom Roberts, all rights reserved.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

http://www.gnu.org/copyleft/gpl.html
*/

#include <memory>
#include <vector>
#include <sstream>

#include "G4VisAttributes.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Color.hh"
#include "G4UserLimits.hh"
#include "G4Material.hh"

#include "BLElement.hh"
#include "BLManager.hh"
#include "BLElementField.hh"
#include "BLGlobalField.hh"
#include "DerivativesSolenoid.hh"
#include "BLKillTrack.hh"

/** Wrap the derivatives solenoid inheriting from the BLElement Field */
class DerivativesSolenoidWrapper : public BLElementField {
public:
	DerivativesSolenoidWrapper() {}
	std::shared_ptr<DerivativesSolenoid> solenoid;
	BLCoordinateTransform global2local;
	G4RotationMatrix* rotation =  nullptr;

	void addFieldValue(const G4double point[4], G4double field[6]) const;
};

void DerivativesSolenoidWrapper::addFieldValue(const G4double point[4], G4double field[6]) const {
	G4double relPoint[4];
	global2local.getLocal(relPoint, point);
	std::vector<double> position = {relPoint[0], relPoint[1], relPoint[2]};
	double time = point[3];
	std::vector<double> bfield(3, 0.0);
	try {
		solenoid->GetFieldValue(position, time, bfield);
	} catch (std::string str) {
		G4Exception("DerivativesSolenoid", str.c_str(), FatalException, "");
	}
	G4ThreeVector g4bfield(bfield[0], bfield[1], bfield[2]);
	if (rotation != nullptr) {
		g4bfield = *rotation * g4bfield;
	}
	field[0] += g4bfield[0];
	field[1] += g4bfield[1];
	field[2] += g4bfield[2];
}

/**	BLCMDsolenoid implements a solenoid.
 *
 *	A solenoid consists of a BLCoil and a current value. It is a
 *	BLElement and can be placed geometrically (multiple times). The
 *	geometrical parameters are in the BLCoil, not the BLCMDsolenoid;
 *	this class uses the variables of its BLCoil to determine the
 *	geometry of the solenoid.
 **/
class BLCMDderivativessolenoid : public BLElement {
public:
    std::shared_ptr<DerivativesSolenoid> solenoid;
	const int maxFourierHarmonic=10;
	double minPhysicalR=-1;
	double maxPhysicalR=-1;
    //FFM parameters
	std::vector<double> harmonicList;
    double cellLength = 0.0;
    //TFM parameters
    double centreLength = 0.0;
    double endLength = 1.0;
    double nominalField = 1.0;
    // material parameters
	G4String material="Vacuum";
	G4String color="1,1,1";
    G4String fieldModelStr;

	/// Default constructor. Defines the command, etc.
	BLCMDderivativessolenoid();

	/// Destructor.
	virtual ~BLCMDderivativessolenoid() { }

	/// Copy constructor.
	BLCMDderivativessolenoid(const BLCMDderivativessolenoid& r);

	/// clone()
	BLElement *clone() { return new BLCMDderivativessolenoid(*this); }

	/// commandName() returns "derivativessolenoid".
	G4String commandName() { return "derivativessolenoid"; }

	/// command() implements the solenoid command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments of the command.
	void defineNamedArgs();
    void defineTFMArgs(); // helper for TanhFieldModel
    void defineFFMArgs(); // helper for FourierFieldModel


	/// construct() - construct the solenoid.
	/// Creates a new SolenoidField and adds it to BLGlobalField.
	void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);
    /// Construct the end field model - Tanh (TFM) or Fourier (FFM) FieldModel
    void constructFieldModel();

	/// getLength() returns the length of the solenoid
	G4double getLength() { return solenoid->length; }

	/// getWidth() returns the outer radius of the solenoid
	G4double getWidth() { return solenoid->maxR*2.0; }

	/// getHeight() returns the outer radius of the solenoid
	G4double getHeight() { return solenoid->maxR*2.0; }

	/// getSurveyPoint() returns points in LOCAL coordinates.
	G4ThreeVector getSurveyPoint(int index) {
		if(index == 0) return G4ThreeVector(0.0, 0.0, -solenoid->length/2.0);
		if(index == 1) return G4ThreeVector(0.0, 0.0, solenoid->length/2.0);
		throw "UNIMPLEMENTED";
	}

	/// isOK() returns true.
	G4bool isOK() { return true; }

	/// generatePoints() from BLElement
	void generatePoints(int npoints, std::vector<G4ThreeVector> &v);

	/// isOutside() from BLElement
	G4bool isOutside(G4ThreeVector &local, G4double tolerance);

};

BLCMDderivativessolenoid defaultDerivativesSolenoid;	// default object


// Default constructor - be sure to use the default constructor BLElement()
BLCMDderivativessolenoid::BLCMDderivativessolenoid() : BLElement(), solenoid(), harmonicList(maxFourierHarmonic)
{
	// register the commandName(), and its synopsis and description.
	registerCommand(BLCMDTYPE_ELEMENT);
	setSynopsis("defines a solenoid (based on an on-axis field and derivatives)");
	setDescription("Solenoid based on field and derivatives.");
	solenoid.reset(new DerivativesSolenoid());
	solenoid->fieldModel.reset(new FourierFieldModel());
}

// Copy constructor - be sure to use the copy constructor BLElement(r)
BLCMDderivativessolenoid::BLCMDderivativessolenoid(const BLCMDderivativessolenoid& r) : BLElement(r), solenoid(r.solenoid->Clone()), harmonicList(r.harmonicList)
{
	// copy fields one at a time (transfers default values from the
	// default object to this new object).
}

int BLCMDderivativessolenoid::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("solenoid: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultDerivativesSolenoid.handleNamedArgs(namedArgs);
	}

	BLCMDderivativessolenoid *t = new BLCMDderivativessolenoid(defaultDerivativesSolenoid);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);
	t->print(argv[0]);

	return retval;
}

void BLCMDderivativessolenoid::defineNamedArgs()
{
	argInt(solenoid->maxDerivative,"maxDerivative","",false);
	argDouble(solenoid->maxR,"maxR","");
	argDouble(solenoid->length,"length","");
	argDouble(minPhysicalR,"minPhysicalR","");
	argDouble(maxPhysicalR,"maxPhysicalR","");
	argString(material, "material", "");
	argString(color, "color", "");
    argString(fieldModelStr, "fieldModel", "");
    argDouble(centreLength,"centreLength","");
    argDouble(endLength,"endLength","");
    argDouble(nominalField, "nominalField", "", 1e-3);
    argDouble(cellLength,"cellLength","");
    defineFFMArgs(); // harmonicList
}

void BLCMDderivativessolenoid::defineFFMArgs() {
    for (size_t i = 0; i < maxFourierHarmonic; ++i) {
        std::stringstream ss;
        ss << "harmonic" << i;
        argDouble(harmonicList[i], ss.str().c_str(), "", 1e-3);
    }
}

void BLCMDderivativessolenoid::constructFieldModel() {
    if (fieldModelStr == "tanh") {
        TanhFieldModel* tfm = new TanhFieldModel();
        tfm->_x0 = centreLength/2;
        tfm->_lambda = endLength;
        tfm->_b0 = nominalField;
        solenoid->fieldModel.reset(tfm);
    } else if (fieldModelStr == "fourier") {
        FourierFieldModel* ffm = new FourierFieldModel();
        ffm->cellLength = cellLength;
        ffm->harmonicList = harmonicList;
        solenoid->fieldModel.reset(ffm);
        std::cerr << "Constructed FFM " << ffm->cellLength << " " << ffm->harmonicList[0] << std::endl;
    } else {
        std::string msg = "Did not recognise fieldModel '"+fieldModelStr+"'";
        G4Exception("BLCMDderivativessolenoid", msg.c_str(), FatalException, "");
    }
}


void BLCMDderivativessolenoid::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)
{
	// geant4 rotation convention is backwards from g4beamline
	G4RotationMatrix *g4rot = 0;
	if(relativeRotation)
		g4rot = new G4RotationMatrix(relativeRotation->inverse());

	// get globalRotation and globalPosition
	G4RotationMatrix *globalRotation = 0;
	if(relativeRotation && parentRotation) {
		globalRotation = 
		    new G4RotationMatrix(*parentRotation * *relativeRotation);
	} else if(relativeRotation) {
		globalRotation = relativeRotation;
	} else if(parentRotation) {
		globalRotation = parentRotation;
	}
	G4ThreeVector globalPosition(relativePosition + parentPosition);
	if(parentRotation)
		globalPosition = *parentRotation * relativePosition +
				parentPosition;

	BLCoordinateTransform global2local(globalRotation, globalPosition);

	G4String thisname = parentName+getName();
	// construct the physical coil
	if (minPhysicalR > 0.0 && maxPhysicalR > 0.0) {
		G4Tubs *tubs = new G4Tubs(thisname+"Tubs", minPhysicalR, 
				maxPhysicalR, getLength()/2.0, 0.0,2.0*pi);
		G4Material *mat = getMaterial(material);;
		G4LogicalVolume *lv = new G4LogicalVolume(tubs,mat, thisname+"LogVol");
		lv->SetVisAttributes(getVisAttrib(color));

		G4PVPlacement *pv = new G4PVPlacement(g4rot,relativePosition,lv,
						thisname, parent,false,0,surfaceCheck);
	}

	DerivativesSolenoidWrapper *ds = new DerivativesSolenoidWrapper();
    try {
        constructFieldModel();
    	solenoid->fieldModel->Initialise(solenoid->maxDerivative);
    } catch (std::string msg) {
        G4Exception("BLCMDderivativessolenoid", msg.c_str(), FatalException, "");
    }
 	ds->solenoid = solenoid;
	ds->global2local = global2local;
	if(global2local.isRotated()) {
		ds->rotation = new G4RotationMatrix(global2local.getRotation());
		ds->rotation->invert();
	}

	BLGlobalField::getObject()->addElementField(ds);

	printf("BLCMDderivativessolenoid::Construct %s parent=%s relZ=%.1f globZ=%.1f\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2],
		global2local.getPosition()[2]);
	printf("\tglobal pos=%.1f,%.1f,%.1f  ",globalPosition[0],
				globalPosition[1],globalPosition[2]);
}

void BLCMDderivativessolenoid::generatePoints(int npoints, std::vector<G4ThreeVector> &v)
{
	v.clear();
	if(minPhysicalR < 0 || maxPhysicalR < 0) return;
	generateTubs(npoints, minPhysicalR, maxPhysicalR, 0.0, 
			360.0*deg, getLength(), v);
}

G4bool BLCMDderivativessolenoid::isOutside(G4ThreeVector &local, G4double tolerance)
{
	if(getLength() == 0) return true;
	G4double r = sqrt(local[0]*local[0]+local[1]*local[1]);
	return r < minPhysicalR+tolerance ||
	       r > maxPhysicalR-tolerance ||
		fabs(local[2]) > getLength()/2.0-tolerance;
}
