package com.sbdm.makeHybridFromHMP;

import java.io.*;
import java.util.*;

/***
 * The tool should be able to process the Single Nucleotide Hapmap file to generte hybrid data.
 * author: s.sivasubramani@cgiar.org
 * date: 12-02-2020
 */
public class Main {

    private static String line;
    private static long lineCount;
    private static ArrayList<String> parents = new ArrayList<>();
    private static ArrayList<String> hybrids = new ArrayList<>();
    private static ArrayList<String> combinedList = new ArrayList<>();
    private static ArrayList<String> calls = new ArrayList<>();
    private static String[] markerLine;
    private static String[] headerLine;
    private static final Map<String, String> iupacHash = new HashMap<String, String>(){{
        put("AA", "A");
        put("TT", "T");
        put("GG", "G");
        put("CC", "C");
        put("AT", "W");
        put("TA", "W");
        put("AG", "R");
        put("GA", "R");
        put("AC", "M");
        put("CA", "M");
        put("TG", "K");
        put("GT", "K");
        put("TC", "Y");
        put("CT", "Y");
        put("GC", "S");
        put("CG", "S");
    }};
    private static final Map<String, String> iupacTwoNuclHash = new HashMap<String, String>(){{
        put("A", "AA");
        put("T", "TT");
        put("G", "GG");
        put("C", "CC");
        put("W", "AT");
        put("R", "AG");
        put("M", "AC");
        put("K", "TG");
        put("Y", "TC");
        put("S", "GC");
        put("N", "NN");
    }};
    private static ArrayList<String> parAlist = new ArrayList<>();
    private static ArrayList<String> parBlist = new ArrayList<>();
    private static int parAcount = 0;
    private static int parBcount = 0;



    public static void main(String[] args) throws IOException {

        if(args.length == 0){
            System.out.println("USAGE: java -jar makeHybridFromHMP.jar input.hmp.txt listOfParentsA.txt listOfParentsB.txt outputPrefix");
            System.exit(0);
        }

        String inHMP = args[0];
        String parA = args[1];
        String parB = args[2];
        String prefix = args[3];
        String outHMP = prefix + ".hmp.txt";
        String outMat = prefix + ".genotypeMatrix.txt";

        /***
         * Reading First Parents list file and creating a List
         */
        try (BufferedReader brParA = new BufferedReader(new FileReader(parA))){
            while ((line = brParA.readLine()) != null){
                parAcount++;
                parAlist.add(line);
            }
            System.out.println("Found " + parAcount + " parent(s) in " + parA + " file.");
        }

        /***
         * Reading Second Parents list file and creating a List
         */
        try (BufferedReader brParB = new BufferedReader(new FileReader(parB))){
            while ((line = brParB.readLine()) != null){
                parBcount++;
                parBlist.add(line);
            }
            System.out.println("Found " + parBcount + " parent(s) in " + parB + " file.");
        }
        combinedList = (ArrayList<String>) getCombinedList(parAlist, parBlist);

        try (BufferedReader brHMP = new BufferedReader(new FileReader(inHMP));
             BufferedWriter bwHMP = new BufferedWriter(new FileWriter(outHMP));
             BufferedWriter bwMat = new BufferedWriter(new FileWriter(outMat))){
            while ((line = brHMP.readLine()) != null){
                lineCount++;
                if((lineCount-1) % 1000 == 0 && lineCount > 1){
                    System.out.println("Processed " + (lineCount-1) + " markers..");
                }
                if(lineCount == 1){
                    headerLine = line.split("\t");
                    for (int i=11; i < headerLine.length; i++){
                        parents.add(headerLine[i]);
                    }
                    for (String parentA: parents) {
                        for (String parentB: parents) {
                            if(!Objects.equals(parentA, parentB) && parAlist.contains(parentA) && parBlist.contains(parentB)){
                                String tempHybridOne = parentA + "x" + parentB;
                                String tempHybridTwo = parentB + "xxxxx" + parentA;
                                if(!(hybrids.contains(tempHybridOne)) && !(hybrids.contains(tempHybridTwo))){
                                    hybrids.add(parentA + "xxxxx" + parentB);
                                }
                            }
                        }
                    }

                    bwHMP.write(getParentsDataHapMap(line, parents, combinedList));
                    bwMat.write("Marker" + "\t" + getParentsDataNumeric(headerLine, parents, combinedList));
                    for (String hybrid:hybrids) {
                        bwHMP.write("\t" + hybrid.replace("xxxxx", "x"));
                        bwMat.write("\t" + hybrid.replace("xxxxx", "x"));
                    }
                    bwHMP.write("\n");
                    bwMat.write("\n");
                }
                else{
//                    String[] alleles = new String[4];
                    markerLine = line.split("\t");
                    ArrayList<String> alleles = getMarkerAlleles(markerLine);
                    if(alleles.size() == 0){
                        System.out.println("Found Missing Data marker. Skipping.. " + markerLine[0] + "\t" + "N");
                        continue;
                    }
                    if(alleles.size() == 1){
                        System.out.println("Found Monomorphic marker. Skipping.. " + markerLine[0] + "\t" + alleles.toString());
                        continue;
                    }
                    if(alleles.size() > 2){
                        System.out.println("Found More than 2 allele containing marker. Skipping.. " + markerLine[0] + "\t" + alleles.toString());
                        continue;
                    }
                    for (int i=11; i < markerLine.length; i++){
                        if(markerLine[i].equals(alleles.get(0))){
                            markerLine[i] = "0";
                        }
                        else if(markerLine[i].equals(alleles.get(1))){
                            markerLine[i] = "2";
                        }
                        else if(markerLine[i].equals("N")){
                            markerLine[i] = "N";
                        }
                        else {
                            markerLine[i] = "1";
                        }
                        calls.add(markerLine[i]);
                    }
                    bwHMP.write(getParentsDataHapMap(line, parents, combinedList));
                    bwMat.write(markerLine[0] + "\t" + getParentsDataNumeric(markerLine, parents, combinedList));
                    for(String hybrid: hybrids){
                        String[] hybridParents = hybrid.split("xxxxx");
                        int indexOfParentA = parents.indexOf(hybridParents[0]);
                        int indexOfParentB = parents.indexOf(hybridParents[1]);
                        String callA = calls.get(indexOfParentA);
                        String callB = calls.get(indexOfParentB);
                        if(Objects.equals(callA, "N") || Objects.equals(callB, "N")){
                            bwHMP.write("\t" + "N");
                            bwMat.write("\t" + "N");
                        }
                        else{
                            float hybridCall = (Float.parseFloat(callA) + Float.parseFloat(callB)) / 2;
                            bwHMP.write("\t" + getAllele(alleles, hybridCall));
                            bwMat.write("\t" + String.valueOf(hybridCall).replace(".0",""));
                        }
                    }
                    bwHMP.write("\n");
                    bwMat.write("\n");
                    calls.clear();
                }
            }
        }
    }

    /***
     * Method to extract Major and Minor alleles from array of HapMap line
     * @param markerLine
     * @return
     */
    private static ArrayList<String> getMarkerAlleles(String[] markerLine) {
        ArrayList<String> alleles = new ArrayList<>();
        String baseOne;
        String baseTwo;
        Map<String, Integer> countAlleles = new HashMap<String, Integer>(){};
        for(int i=11;i<markerLine.length;i++){
            if(!iupacTwoNuclHash.containsKey(markerLine[i])){
                System.out.println("Invalid or Polyploid data \"" + markerLine[i] + "\" found for Marker " + markerLine[0] + ". Please check the data.");
            }
            baseOne = String.valueOf(iupacTwoNuclHash.get(markerLine[i]).charAt(0));
            baseTwo = String.valueOf(iupacTwoNuclHash.get(markerLine[i]).charAt(1));
            if(countAlleles.containsKey(baseOne)){
                countAlleles.put(baseOne, countAlleles.get(baseOne)+1);
            }
            else{
                countAlleles.put(baseOne, 1);
            }
            if(countAlleles.containsKey(baseTwo)){
                countAlleles.put(baseTwo, countAlleles.get(baseTwo)+1);
            }
            else{
                countAlleles.put(baseTwo, 1);
            }
        }
        LinkedHashMap<String, Integer> reverseSortedMap = new LinkedHashMap<>();
        countAlleles.entrySet()
                .stream()
                .sorted(Map.Entry.comparingByValue(Comparator.reverseOrder()))
                .forEachOrdered(x -> reverseSortedMap.put(x.getKey(), x.getValue()));
        int index = 0;
        for (Map.Entry<String, Integer> count : reverseSortedMap.entrySet()) {
            if(!count.getKey().equals("N")){
                alleles.add(count.getKey());
                index++;
            }
        }
        if(alleles.size() > 1) {
            if (countAlleles.get(alleles.get(0)).equals(countAlleles.get(alleles.get(1)))) {
                Collections.sort(alleles);
            }
        }
        return alleles;
    }

    /***
     * Method to return IUPAC call from numberic
     * @param alleles
     * @param hybridCall
     * @return
     */
    private static String getAllele(ArrayList<String> alleles, float hybridCall) {
        String allele;
        if(hybridCall == 0){
            allele = alleles.get(0);
        }
        else if(hybridCall == 2){
            allele = alleles.get(1);
        }
        else{
            allele = iupacHash.get(alleles.get(0) + alleles.get(1));
        }
        return allele;
    }

    /***
     * Method to return Numeric data for Parents
     * @param arrayLine
     * @param parents
     * @param combinedList
     * @return
     */
    private static String getParentsDataNumeric(String[] arrayLine, ArrayList<String> parents, ArrayList<String> combinedList) {
        StringBuilder parentsLine = new StringBuilder();
        for (String parent : combinedList) {
            if(parentsLine.toString().equals("")){
                parentsLine = new StringBuilder(arrayLine[parents.indexOf(parent) + 11]);
            }
            else {
                parentsLine.append("\t").append(arrayLine[parents.indexOf(parent) + 11]);
            }
        }
        return parentsLine.toString();
    }

    /***
     * Method to return HapMap data for parents
     * @param line
     * @param parents
     * @param combinedList
     * @return
     */
    private static String getParentsDataHapMap(String line, ArrayList<String> parents, ArrayList<String> combinedList) {
        StringBuilder parentsLine = new StringBuilder();
        String[] arrayLine = line.split("\t");
        for (int i=0; i < 11; i++){
            if(parentsLine.toString().equals("")) {
                parentsLine = new StringBuilder(arrayLine[i]);
            }
            else {
                parentsLine.append("\t").append(arrayLine[i]);
            }
        }
        for (String parent : combinedList) {
            parentsLine.append("\t").append(arrayLine[parents.indexOf(parent) + 11]);
        }
        return parentsLine.toString();
    }

    /***
     * Method to return combined unique list after merging 2 lists
     * @param parAlist
     * @param parBlist
     * @return
     */
    private static List<String> getCombinedList(ArrayList<String> parAlist, ArrayList<String> parBlist) {
        Set<String> set = new LinkedHashSet<>(parAlist);
        set.addAll(parBlist);
        return new ArrayList<>(set);
    }
}
