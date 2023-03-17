package fitnessruch2bndlng;

/*
 *   Copyright (C) 2016
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Affero General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Affero General Public License for more details.
 *
 *   You should have received a copy of the GNU Affero General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.vecmath.Point3d;

import org.openscience.cdk.ChemFile;
import org.openscience.cdk.ChemObject;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.smsd.Isomorphism;
import org.openscience.cdk.smsd.interfaces.Algorithm;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;

import denoptim.io.DenoptimIO;
import denoptim.utils.DENOPTIMMathUtils;
import denoptim.utils.DoubleVector;
import denoptim.utils.GenUtils;

/**
 * Tool for the analysis of Ru 14-electrons compounds as catalysts for 
 * olefine metathesis. 
 * <br><br>
 * WARNING! Only complexes reflecting the formula Ru(Cl)(Cl)(L)=CH2
 * will be evaluated by this tool.
 * <br><br>
 * FITNESS=-(Ru=CH2 bond length)
 * <br><br>
 * The task includes steps as (i) evaluation of geometry constraints,
 * and (ii) calculation of the fitness
 * <br><br>
 * 
 * @author Vishwesh Venkatraman
 * @author Marco Foscato
 * @author Sondre Eliasson
 */

public class FitnessRuCH2BndLng
{

    private static final Logger LOGGER = Logger.getLogger(
            FitnessRuCH2BndLng.class.getName());

    //SDF pre-geometry optimization with MOPAC
    private static String inpSdfFile;

    //SDF optimized geometry
    private static String optSdfFile;

    //XYZ file (coordinates with higher precision than SDF)
    private static String hpXYZ;

    //SDF output file to be produced
    private static String outsdfFile;

    private static double MinNonBondedDistance = Double.MAX_VALUE;
    private static double MaxBondedDistance = Double.MIN_VALUE;
    private static double MaxTorsion = Double.MIN_VALUE;
    private static double MinAngle = Double.MAX_VALUE;
    private static String wrkDir;

//------------------------------------------------------------------------------    

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args)
    {
        if (args.length < 1)
        {
            System.err.println("Usage: java -jar FitnessRuCH2BndLng.jar parameterFile");
            System.exit(-1);
        }

        String paramFile = args[0];

        try
        {
            readParameters(paramFile);
            checkParameters();

            // Read graph representation (so-called 2D)            
            IAtomContainer mol = DenoptimIO.readSingleSDFFile(inpSdfFile);

            // Read Gaussian optimized geometry
            IAtomContainer optMol = DenoptimIO.readSingleSDFFile(optSdfFile);

            // Read higher precision XYZ coordinates
            ArrayList<String> xyzTxt = DenoptimIO.readList(hpXYZ);
            if (xyzTxt.size() != (optMol.getAtomCount() + 2))
            {
                    LOGGER.log(Level.SEVERE, "Inconsistency between SDF and XYZ: check output from DFT");
                    System.exit(-1);
            }
            for (int il=2; il<xyzTxt.size(); il++)
            {
            String line = xyzTxt.get(il);
            line = line.trim();
            String[] parts = line.split("\\s+");
    
            String sym = parts[0];
            double hpx = Double.parseDouble(parts[1]); 
                    double hpy = Double.parseDouble(parts[2]);
                    double hpz = Double.parseDouble(parts[3]);
            Point3d hp3d = new Point3d(hpx,hpy,hpz);
    
            IAtom a = optMol.getAtom(il-2);
    
            //Check that the atom is the same
            boolean isok = true;
            if (!a.getSymbol().equals(sym))
            {
                LOGGER.log(Level.SEVERE, "Inconsistency between SDF and XYZ: check atom " + a);
                        System.exit(-1); 
            }
                    if (hp3d.distance(a.getPoint3d()) > 0.0003)
            {
                        LOGGER.log(Level.SEVERE, "Inconsistency between SDF and XYZ: check position of atom " + a);
                        System.exit(-1);           
            }
    
            //read in higher precision coordinates
            a.setPoint3d(hp3d);
            }

            // Copy SDF properties to optMol
            Map<Object,Object> molProps = mol.getProperties();
            for (Object key : molProps.keySet())
            {
                optMol.setProperty(key,molProps.get(key));
            }

            // Define indexes of atoms in Ru(L)(Cl)(Cl)=CH2
            Map<String,Integer> atomIndeces = defineAtomIndexes(optMol);
        
            //Calculate descriptors
            DoubleVector descriptors = new DoubleVector(7);
            calculateDescriptors(optMol, atomIndeces, descriptors);
            StringBuilder sb = new StringBuilder();
            for (int j=0; j<descriptors.length(); j++)
            {
                sb.append(String.format("%6.5f", 
                                         descriptors.getValue(j))).append(" ");
            }
            optMol.setProperty("Descriptors", sb.toString().trim());

            // Check constraints
             String status = checkMolecule(descriptors, atomIndeces, optMol);
            if (!status.equalsIgnoreCase("OK"))
            {
                // write MOL_ERROR tag
                String msg = "#threshold violated: " + status;
                optMol.setProperty("MOL_ERROR", status);
                DenoptimIO.writeMolecule(outsdfFile, optMol, false);
            }
            else
            {
                optMol.setProperty("FITNESS",String.format("%8.5f", 
                              GenUtils.roundValue(-1.0 * descriptors.getValue(4), 5)));
                optMol.setProperty("calculated_ATOM_INDECES",atomIndeces);
                DenoptimIO.writeMolecule(outsdfFile, optMol, false);
            }
        }
        catch (Exception de)
        {
            LOGGER.log(Level.SEVERE, null, de);
            System.exit(-1);
        }

        System.exit(0);
    }

//------------------------------------------------------------------------------
    /**
     * Looks for the indexes of the atoms involved in the calculation of 
     * descriptors and the constraints.
     *
     * WARNING! This assumes Ru compounds with Ru(Cl)(Cl)(L)=CH2 core
     */
    private static Map<String,Integer> defineAtomIndexes(IAtomContainer mol)
                throws Exception
    {
        Map<String,Integer> atmIndeces = new HashMap<String,Integer>();

        //Ru
        int iRu = -1;
        for (IAtom atm : mol.atoms())
        {
            if (atm.getSymbol().equals("Ru"))
            {
                iRu = mol.getAtomNumber(atm);
            }
        }
        if (iRu < 0)
        {
            LOGGER.log(Level.SEVERE, "Failed to identify Ru atom. Check file " + 
                                                                         outsdfFile);
            DenoptimIO.writeMolecule(outsdfFile, mol, false);
            System.exit(-1);
        }

        atmIndeces.put("iRu",iRu);
        atmIndeces.put("iCl1",-1);
        atmIndeces.put("iCl2",-1);
        atmIndeces.put("iC",-1);
        atmIndeces.put("iL",-1);
        atmIndeces.put("iH1",-1);
        atmIndeces.put("iH2",-1);
        List<IAtom> nbrsRu = mol.getConnectedAtomsList(mol.getAtom(iRu));
        for (IAtom nbr : nbrsRu)
        {
            //Cl
            if (nbr.getSymbol().equals("Cl"))
            {
                if (atmIndeces.get("iCl1") < 0)
                {
                    atmIndeces.put("iCl1",mol.getAtomNumber(nbr));
                } else {
                    atmIndeces.put("iCl2",mol.getAtomNumber(nbr));
                }
                continue;
            } 

            //C and H in methylidene
            if (nbr.getSymbol().equals("C"))
            {
                List<IAtom> nbrsC = mol.getConnectedAtomsList(nbr);
                if (nbrsC.size() == 3)
                {                
                    int tmp1 = -1;
                    for (IAtom nbrC : nbrsC)
                    {
                        if (nbrC.getSymbol().equals("H"))
                        {
                            if (tmp1 < 0)
                            {
                                tmp1 = mol.getAtomNumber(nbrC);
                            } else {
                                atmIndeces.put("iC",mol.getAtomNumber(nbr));
                                atmIndeces.put("iH1",tmp1);
                                atmIndeces.put("iH2",mol.getAtomNumber(nbrC));
                            }
                        }
                    } // end loop over nbrs of C
                }
            }
        } //end loop over nbrs of Ru

        //L-ligand
        for (IAtom nbr : nbrsRu)
        {
            int tmp2 = mol.getAtomNumber(nbr);

            if (tmp2 == atmIndeces.get("iC"))
                continue;

            if (tmp2 == atmIndeces.get("iCl1"))
                continue;

            if (tmp2 == atmIndeces.get("iCl2"))
                continue;

            atmIndeces.put("iL",tmp2);
        }

        //Check indeces
        for (String label : atmIndeces.keySet())        
        {
            //Check assigniation
            if (atmIndeces.get(label) < 0)
            {
                LOGGER.log(Level.SEVERE, "Failed to identify '" + label
                                + "' atom. Check file " + outsdfFile);
                DenoptimIO.writeMolecule(outsdfFile, mol, false);
                System.exit(-1);                
            }
        
            //Check unicity
            for (String label2 : atmIndeces.keySet())
            {
        if (label.equals(label2))
            continue;

                if (atmIndeces.get(label).equals(atmIndeces.get(label2)))
                {
                    LOGGER.log(Level.SEVERE, "Ambiguity in identifying '" + label
                                + "' and '" + label2 + "'. Check file " + outsdfFile);
                    DenoptimIO.writeMolecule(outsdfFile, mol, false);
                    System.exit(-1);
                }
            }
        }

        return atmIndeces;
    }
    
//------------------------------------------------------------------------------
    
    private static void checkParameters() throws Exception
    {
        if (inpSdfFile == null || inpSdfFile.length() == 0)
        {
            throw new Exception("Input SDF file not supplied. Check parameter file.");
        }
        if (hpXYZ == null || hpXYZ.length() == 0)
        {
            throw new Exception("Input XYZ  file not supplied. Check parameter file.");
        }
        if (wrkDir == null || wrkDir.length() == 0)
        {
            throw new Exception("Working directory not supplied. Check parameter file.");
        }
        if (outsdfFile == null || outsdfFile.length() == 0)
        {
            throw new Exception("Output SDF file not supplied. Check parameter file.");
        }
        if (optSdfFile == null || optSdfFile.length() == 0)
        {
            throw new Exception("Optimized geometry file not supplied. Check parameter file.");
        }
    }

//------------------------------------------------------------------------------

    /**
     * Calculate descriptors for the molecule
     *
     * @param mol
     * @param atomIndeces map with the indeces of the atoms to be used
     * @param descriptors
     * @param atom_charges
     * @throws Exception
     */
    private static void calculateDescriptors(IAtomContainer mol, 
            Map<String,Integer> atomIndeces,
            DoubleVector descriptors)
            throws Exception
    {
        int h1 = atomIndeces.get("iH1");
        int h2 = atomIndeces.get("iH2");
        int Cl1 = atomIndeces.get("iCl1");
        int Cl2 = atomIndeces.get("iCl2");
        int Ru = atomIndeces.get("iRu");
        int C = atomIndeces.get("iC");
        int X = atomIndeces.get("iL");

        double distRuCl = getDistance(mol, Ru, Cl1);
        distRuCl += getDistance(mol, Ru, Cl2);
        distRuCl /= 2;
        descriptors.setValue(0, distRuCl);

        double angleClRuCl = getAngle(mol, Cl1, Ru, Cl2);
        descriptors.setValue(1, angleClRuCl);

        double angleClRuC = getAngle(mol, Cl1, Ru, C);
        angleClRuC += getAngle(mol, Cl2, Ru, C);
        angleClRuC /= 2;
        descriptors.setValue(2, angleClRuC);

        double d1 = Math.abs(getTorsion(mol, h1, C, Ru, X));
        double d2 = Math.abs(getTorsion(mol, h2, C, Ru, X));
        double torsionXRuCH = Math.min(Math.abs(d1), Math.abs(d2));
        descriptors.setValue(3, torsionXRuCH);

        double distRuC = getDistance(mol, Ru, C);
        descriptors.setValue(4, distRuC);

        double distRuL = getDistance(mol, Ru, X);
        descriptors.setValue(5, distRuL);

        double angleCRuL = getAngle(mol, C, Ru, X);
        descriptors.setValue(6, angleCRuL);

    }

//------------------------------------------------------------------------------

    /**
     * Calculate the minimum non-bonded Distance
     *
     * @param scaffold
     * @param mol
     * @return the distance
     */
    private static double getMinimumNonBondedAtomDistance(
                        Map<String,Integer> atomIndeces,
                            IAtomContainer mol)
    {
        int Ru = atomIndeces.get("iRu");

        IAtom atmRu = mol.getAtom(Ru);

        //List of atoms directly connected to Ru 
        ArrayList<IAtom> lst = new ArrayList<>();
        List<IAtom> nbrsRu = mol.getConnectedAtomsList(atmRu);
        lst.addAll(nbrsRu);

        //Add also all the atoms involved in calculation of descriptors
        for (String label : atomIndeces.keySet())
        {
            IAtom atm = mol.getAtom(atomIndeces.get(label));
            if (!lst.contains(atm))
            {
                lst.add(atm);
            }
        }

        // measure distances for all non-bonded atoms
        double d = Double.MAX_VALUE;

        int k = 0;
        for (IAtom atom : mol.atoms())
        {
            if (!lst.contains(atom))
            {
                double f = getDistance(mol, Ru, k);
                if (f < d)
                {
                    d = f;
                }
            }
            k++;
        }
        return d;
    }

//------------------------------------------------------------------------------

    private static double getDistance(IAtomContainer mol, int i, int j)
    {
        Point3d p1 = mol.getAtom(i).getPoint3d();
        Point3d p2 = mol.getAtom(j).getPoint3d();
        return DENOPTIMMathUtils.distance(p1, p2);
    }

//------------------------------------------------------------------------------

    private static double getAngle(IAtomContainer mol, int i, int j, int k)
    {
        Point3d p1 = mol.getAtom(i).getPoint3d();
        Point3d p2 = mol.getAtom(j).getPoint3d();
        Point3d p3 = mol.getAtom(k).getPoint3d();
        return DENOPTIMMathUtils.angle(p1, p2, p3);
    }

//------------------------------------------------------------------------------

    private static double getTorsion(IAtomContainer mol, int i, int j, int k, int l)
    {
        Point3d p1 = mol.getAtom(i).getPoint3d();
        Point3d p2 = mol.getAtom(j).getPoint3d();
        Point3d p3 = mol.getAtom(k).getPoint3d();
        Point3d p4 = mol.getAtom(l).getPoint3d();
        return DENOPTIMMathUtils.torsion(p1, p2, p3, p4);
    }

//------------------------------------------------------------------------------

    /**
     *
     * @param descriptors
     * @param scaffold
     * @param mol
     * @return the error message due to constraint violation or "OK" if the
     * constraints are satisfied
     */
    private static String checkMolecule(DoubleVector descriptors,
                               Map<String,Integer> atomIndeces,
                               IAtomContainer mol) throws Exception
    {
        //threshold_dist=Ru-X|(0,2.5); X-X|(0,2.5); Ru~X@n|(2.9,)
        //threshold_angle=Cl-Ru-C|(90,)
        //threshold_dihedral=H-C-Ru-X@n|(,10)
        String msg = null;
        if (descriptors.getValue(0) > MaxBondedDistance)
        {
            msg = ": Maximum bonded distance threshold violated. "
                    + String.format("%4f", descriptors.getValue(0));
            LOGGER.info(msg);
            return "Maximum bonded distance " + 
                    String.format("%2f", descriptors.getValue(0));
        }

        if (descriptors.getValue(2) < MinAngle)
        {
            msg = ": Minimum angle threshold violated. "
                    + String.format("%4f", descriptors.getValue(2));
            LOGGER.info(msg);
            return "Minimum angle " + 
                    String.format("%2f", descriptors.getValue(2));
        }

        if (descriptors.getValue(3) > MaxTorsion)
        {
            msg = ": Torsion angle threshold violated. "
                    + String.format("%4f", descriptors.getValue(3));
            LOGGER.info(msg);
            return "Torsion angle " + 
                    String.format("%2f", descriptors.getValue(3));
        }

        double f = getMinimumNonBondedAtomDistance(atomIndeces,mol);
        if (f < MinNonBondedDistance)
        {
            msg = ": Minimum non-bonded distance threshold violated. "
                    + String.format("%4f", f);
            LOGGER.info(msg);
            return "Minimum non-bonded distance " + 
                    String.format("%2f", f);
        }

        return "OK";
    }

//------------------------------------------------------------------------------

    private static void readParameters(String filename) throws Exception
    {
        BufferedReader br = null;
        String option, line;

        try
        {
            br = new BufferedReader(new FileReader(filename));
            while ((line = br.readLine()) != null)
            {
                line = line.trim();
                if (line.length() == 0)
                {
                    continue;
                }

                if (line.startsWith("#"))
                {
                    continue;
                }

                if (line.toUpperCase().startsWith("INPSDF"))
                {
                    option = line.substring(line.indexOf("=") + 1).trim();
                    if (option.length() > 0)
                    {
                        inpSdfFile = option;
                    }
                    continue;
                }

                if (line.toUpperCase().startsWith("OPTSDF"))
                {
                    option = line.substring(line.indexOf("=") + 1).trim();
                    if (option.length() > 0)
                    {
                        optSdfFile = option;
                    }
                    continue;
                }
                
                if (line.toUpperCase().startsWith("OUTSDF"))
                {
                    option = line.substring(line.indexOf("=") + 1).trim();
                    if (option.length() > 0)
                    {
                        outsdfFile = option;
                    }
                    continue;
                }

                if (line.toUpperCase().startsWith("MAXBNDDIST"))
                {
                    option = line.substring(line.indexOf("=") + 1).trim();
                    if (option.length() > 0)
                    {
                        MaxBondedDistance = Double.parseDouble(option);
                    }
                    continue;
                }

                if (line.toUpperCase().startsWith("MINANGLE"))
                {
                    option = line.substring(line.indexOf("=") + 1).trim();
                    if (option.length() > 0)
                    {
                        MinAngle = Double.parseDouble(option);
                    }
                    continue;
                }

                if (line.toUpperCase().startsWith("MAXTORSION"))
                {
                    option = line.substring(line.indexOf("=") + 1).trim();
                    if (option.length() > 0)
                    {
                        MaxTorsion = Double.parseDouble(option);
                    }
                    continue;
                }

                if (line.toUpperCase().startsWith("MINNBDIST"))
                {
                    option = line.substring(line.indexOf("=") + 1).trim();
                    if (option.length() > 0)
                    {
                        MinNonBondedDistance = Double.parseDouble(option);
                    }
                    continue;
                }

                if (line.toUpperCase().startsWith("WORKDIR"))
                {
                    option = line.substring(line.indexOf("=") + 1).trim();
                    if (option.length() > 0)
                    {
                        wrkDir = option;
                    }
                    continue;
                }

                if (line.toUpperCase().startsWith("HPXYZ"))
                {
                    option = line.substring(line.indexOf("=") + 1).trim();
                    if (option.length() > 0)
                    {
                        hpXYZ = option;
                    }
                    continue;
                }
            }
        }
        finally
        {
            if (br != null)
            {
                br.close();
            }
        }
    }

//------------------------------------------------------------------------------

    private static IAtomContainer readSDF(String filename)
            throws Exception
    {
        MDLV2000Reader mdlreader = null;
        ArrayList<IAtomContainer> lstContainers = new ArrayList<>();

        try
        {
            mdlreader = new MDLV2000Reader(new FileReader(new File(filename)));
            ChemFile chemFile = (ChemFile) mdlreader.read((ChemObject) new ChemFile());
            lstContainers.addAll(
                    ChemFileManipulator.getAllAtomContainers(chemFile));
        }
        finally
        {
            if (mdlreader != null)
            {
                mdlreader.close();
            }
        }

        if (lstContainers.isEmpty())
        {
            throw new Exception("No data found in " + filename);
        }

        return lstContainers.get(0);
        //return lstContainers.toArray(new IAtomContainer[lstContainers.size()]);
    }
    
    
//------------------------------------------------------------------------------    
   
}
