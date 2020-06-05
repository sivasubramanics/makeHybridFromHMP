package com.sbdm.cleanReads;

import org.apache.commons.cli.*;

import java.io.*;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.zip.GZIPInputStream;

public class Main {


    private static int noReads = 0;
    private static String readStatus;

    public static void main(String[] args) throws FileNotFoundException {
        // write your code here
        Options options = Utils.getOptions();

        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd = null;

        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.out.println(e.getMessage());
            formatter.printHelp("java -jar clearFastq.jar -f read1.fq -r reaf2.fq -p output",
                    options);
            System.exit(1);
        }


        String inFastqOne = cmd.getOptionValue("forward");
        String inFastqTwo = cmd.getOptionValue("reverse");
        String outPrefix = cmd.getOptionValue("prefix");
        String outFastqOne = outPrefix + ".cleaned.1.fastq";
        String outFastqTwo = outPrefix + ".cleaned.2.fastq";
        String outFastqOrp = outPrefix + ".cleaned.orphan.fastq";
        String outFastqRej = outPrefix + ".cleaned.rejected.fastq";
//        int nThreads = Integer.parseInt(cmd.getOptionValue("threads", "1"));

        try(BufferedReader brOne = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inFastqOne))));
            BufferedReader brTwo = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inFastqTwo))));
            BufferedWriter bwOne = new BufferedWriter(new FileWriter(outFastqOne));
             BufferedWriter bwTwo = new BufferedWriter(new FileWriter(outFastqTwo));
             BufferedWriter bwOrp = new BufferedWriter(new FileWriter(outFastqOrp));
             BufferedWriter bwRej = new BufferedWriter(new FileWriter(outFastqRej));){

            long start;
            long end;
            long duration = 0;

            start = System.nanoTime();
            boolean fqLeft = true;
            while (fqLeft){
                ThreadPoolExecutor executor = (ThreadPoolExecutor) Executors.newFixedThreadPool(2);

                HashMap<Integer, List<String>> fqOneMap = new HashMap<>();
                HashMap<Integer, List<String>> fqTwoMap = new HashMap<>();
                HashMap<Integer, List<String>> fqOneTempMap = new HashMap<>();
                HashMap<Integer, List<String>> fqTwoTempMap = new HashMap<>();
                HashMap<Integer, String> pairStatus = new HashMap<>();


//                List <String> fqOne  = new ArrayList<>();;
//                List <String> fqTwo  = new ArrayList<>();;
//                List <String> fqOneTemp = new ArrayList<String>(fqOne);
//                List <String> fqTwoTemp  = new ArrayList<String>(fqTwo);;

                for (int i = 0; i <= 9; i++){
                    fqOneMap.put(i, Utils.getLine(brOne));
                    fqTwoMap.put(i, Utils.getLine(brTwo));
                }

//                fqOne = Utils.getLine(brOne);
//                fqTwo = Utils.getLine(brTwo);
                Utils.copyMap(fqOneTempMap, fqOneMap);
                Utils.copyMap(fqTwoTempMap, fqTwoMap);

                for(int i=0; i<= 9; i++){
                    List<String> fqOne = fqOneMap.get(i);
                    List<String> fqTwo = fqTwoMap.get(i);
                    List<String> fqOneTemp = fqOneTempMap.get(i);
                    List<String> fqTwoTemp = fqTwoTempMap.get(i);
                    fqLeft = isFqLeft(fqOne, fqTwo, duration, fqLeft);
                    if(!fqLeft)
                        break;
                    if(fqOne.get(1).length() != fqOne.get(3).length()){
                        System.out.println("Irregular Fastq file.\n");
                        System.exit(1);
                    }
                    if(fqTwo.get(1).length() != fqTwo.get(3).length()){
                        System.out.println("Irregular Fastq file.\n");
                        System.exit(1);
                    }

                    pairStatus.put(i, new ProcessFastqRecord(bwOne, bwTwo, bwOrp, bwRej, fqOne, fqTwo, fqOneTemp, fqTwoTemp).invoke());
                    noReads++;
                }


                end = System.nanoTime();
                duration = (end - start) / 1000000000;
                System.out.print("Processed " + noReads + " read pairs... " + duration + "secs" + "\r");
                System.out.flush();
            }

        } catch (IOException e) {
            e.printStackTrace();
        } ;
    }

    public class ProcessThread implements Runnable {
        private ProcessFastqRecord processFastqRecord;

        ProcessThread(ProcessFastqRecord processFastqRecord){
            this.processFastqRecord =processFastqRecord;
        }

        @Override
        public void run(){
            try {
                processFastqRecord.invoke();
            } catch (Exception e) {
                System.out.println("ReaderThread: " + "Error processing file read");
            }
        }
    }

    private static boolean isFqLeft(List<String> fqOne, List<String> fqTwo, long duration, boolean fqLeft) {
        if(fqOne.get(0) == null || fqTwo.get(0) == null){
            System.out.println("Total reads processed: " + noReads + " (" + duration + "secs)");
            fqLeft = false;
        }
        return fqLeft;
    }

    private static String processFastqPair(BufferedWriter bwOne, BufferedWriter bwTwo, BufferedWriter bwOrp, BufferedWriter bwRej, List<String> fqOne, List<String> fqTwo, List<String> fqOneTemp, List<String> fqTwoTemp) throws IOException {
        fqOne.set(3, CoreFuntions.doBtrim(fqOne.get(3)));
        fqTwo.set(3, CoreFuntions.doBtrim(fqTwo.get(3)));

        int startBaseOne = CoreFuntions.doNtrim(fqOne.get(0));
        int startBaseTwo = CoreFuntions.doNtrim(fqTwo.get(0));

        fqOne.set(3, Utils.getHQseq(fqOne.get(3)));
        fqTwo.set(3, Utils.getHQseq(fqTwo.get(3)));

        fqOne.set(1, CoreFuntions.doAdaptCheck(fqOne.get(1)));
        fqTwo.set(1, CoreFuntions.doAdaptCheck(fqTwo.get(1)));

        if(startBaseOne > startBaseTwo){
            fqOne = Utils.updateFqEntry(fqOne, startBaseOne);
            fqTwo = Utils.updateFqEntry(fqTwo, startBaseOne);
        }
        else {
            fqOne = Utils.updateFqEntry(fqOne, startBaseTwo);
            fqTwo = Utils.updateFqEntry(fqTwo, startBaseTwo);
        }

        readStatus = Utils.writeFqEntry(fqOne, fqTwo, fqOneTemp, fqTwoTemp, bwOne, bwTwo, bwOrp, bwRej);
        return readStatus;
    }

    private static class ProcessFastqRecord {
        private BufferedWriter bwOne;
        private BufferedWriter bwTwo;
        private BufferedWriter bwOrp;
        private BufferedWriter bwRej;
        private List<String> fqOne;
        private List<String> fqTwo;
        private List<String> fqOneTemp;
        private List<String> fqTwoTemp;

        public ProcessFastqRecord(BufferedWriter bwOne, BufferedWriter bwTwo, BufferedWriter bwOrp, BufferedWriter bwRej, List<String> fqOne, List<String> fqTwo, List<String> fqOneTemp, List<String> fqTwoTemp) {
            this.bwOne = bwOne;
            this.bwTwo = bwTwo;
            this.bwOrp = bwOrp;
            this.bwRej = bwRej;
            this.fqOne = fqOne;
            this.fqTwo = fqTwo;
            this.fqOneTemp = fqOneTemp;
            this.fqTwoTemp = fqTwoTemp;
        }

        public String invoke() throws IOException {
            return processFastqPair(bwOne, bwTwo, bwOrp, bwRej, fqOne, fqTwo, fqOneTemp, fqTwoTemp);
        }
    }
}
