package com.sbdm.makeHybridFromHMP;

import java.io.*;
import java.util.*;

public class Main {

    public static String line;
    public static long lineCount;
    public static ArrayList<String> parents = new ArrayList<String>();
    public static ArrayList<String> hybrids = new ArrayList<String>();
    public static ArrayList<String> combinedList = new ArrayList<String>();
    public static ArrayList<String> calls = new ArrayList<String>();
    public static String[] markerLine;
    public static String[] headerLine;
    public static Map<String, String> iupacHash = new HashMap<String, String>(){{
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
    public static Map<String, String> iupacTwoNuclHash = new HashMap<String, String>(){{
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
    public static ArrayList<String> parAlist = new ArrayList<>();
    public static ArrayList<String> parBlist = new ArrayList<>();
    public static int parAcount = 0;
    public static int parBcount = 0;



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

        try (BufferedReader brParA = new BufferedReader(new FileReader(parA))){
            while ((line = brParA.readLine()) != null){
                parAcount++;
                parAlist.add(line);
            }
            brParA.close();
        }

        try (BufferedReader brParB = new BufferedReader(new FileReader(parB))){
            while ((line = brParB.readLine()) != null){
                parBcount++;
                parBlist.add(line);
            }
            brParB.close();
        }


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
                            if(parentA != parentB && parAlist.contains(parentA) && parBlist.contains(parentB)){
                                String tempHybridOne = parentA + "x" + parentB;
                                String tempHybridTwo = parentB + "xxxxx" + parentA;
                                if(!(hybrids.contains(tempHybridOne)) && !(hybrids.contains(tempHybridTwo))){
                                    hybrids.add(parentA + "xxxxx" + parentB);
                                }
                            }
                        }
                    }
                    combinedList = (ArrayList<String>) getCombinedList(parAlist, parBlist);
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
                        if(callA == "N" || callB == "N"){
                            bwHMP.write("\t" + "N");
                            bwMat.write("\t" + "N");
                        }
                        else{
                            float hybridCall = (Float.valueOf(callA) + Float.valueOf(callB)) / 2;
                            bwHMP.write("\t" + getAllele(alleles, hybridCall));
                            bwMat.write("\t" + String.valueOf(hybridCall).replace(".0",""));
                        }
                    }
                    bwHMP.write("\n");
                    bwMat.write("\n");
                    calls.clear();
                }
            }
            bwHMP.close();
            bwMat.close();
        }
    }


    private static ArrayList<String> getMarkerAlleles(String[] markerLine) {
        ArrayList<String> alleles = new ArrayList<>();
        String baseOne = "";
        String baseTwo = "";
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

    private static String getAllele(ArrayList<String> alleles, float hybridCall) {
        String allele = "";
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

    private static String getParentsDataNumeric(String[] arrayLine, ArrayList<String> parents, ArrayList<String> combinedList) {
        String parentsLine = "";
        for (String parent : combinedList) {
            if(parentsLine.equals("")){
                parentsLine = arrayLine[parents.indexOf(parent) + 11];
            }
            else {
                parentsLine = parentsLine + "\t" + arrayLine[parents.indexOf(parent) + 11];
            }
        }
        return parentsLine;
    }

    private static String getParentsDataHapMap(String line, ArrayList<String> parents, ArrayList<String> combinedList) {
        String parentsLine = "";
        String[] arrayLine = line.split("\t");
        for (int i=0; i < 11; i++){
            if(parentsLine.equals("")) {
                parentsLine = arrayLine[i];
            }
            else {
                parentsLine = parentsLine + "\t" + arrayLine[i];
            }
        }
        for (String parent : combinedList) {
            parentsLine = parentsLine + "\t" + arrayLine[parents.indexOf(parent) + 11];
        }
        return parentsLine;
    }

    private static List<String> getCombinedList(ArrayList<String> parAlist, ArrayList<String> parBlist) {
        Set<String> set = new LinkedHashSet<>(parAlist);
        set.addAll(parBlist);
        return new ArrayList<>(set);
    }

}
