/*
 * Copyright (C) 2015 BioDAG - Athanassios Kintsakis
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package impute2chromopainter;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;

//import org.apache.commons.io.FileUtils;
/**
 *
 * @author thanos
 */
public class Impute2chromopainter {

    /**
     * @param args the command line arguments
     */
    public static void withBitSet(String filename, String outfile, String geneticMap, String recombOut) throws FileNotFoundException, IOException {
        BufferedReader r1 = new BufferedReader(new FileReader(new File(filename)));
        BufferedReader r2 = new BufferedReader(new FileReader(new File(geneticMap)));
        BufferedWriter wr2 = new BufferedWriter(new FileWriter(new File(recombOut)));
        wr2.write("start.pos recom.rate.perbp");
        wr2.newLine();

        String line;
        String[] tmp;
        int count = 0;
        ArrayList<Integer> pValues = new ArrayList<Integer>();
        ArrayList<BitSet> sets = new ArrayList<BitSet>();
        int lineLength = 0;
        pValues.add(-1);        

        //StringBuilder sb = new StringBuilder();
        while ((line = r1.readLine()) != null) {
            tmp = line.split(" ");
            lineLength = tmp.length;
            int tempConvert = Integer.valueOf(tmp[2]);
            if (tempConvert > (pValues.get(count))) {
                pValues.add(tempConvert);

                BitSet tk = new BitSet();
                for (int i = 5; i < lineLength; i++) {
                    if (tmp[i].equals("1")) {
                        tk.set(i - 5);
                    }
                }
                sets.add(tk);
                count++;
            }
//            else
//            {
//                System.out.println("WARNING: genetic map error -- position duplicated in haplotype input file "+count+" = "+(count+1)+" . Will use only first instance of position.");
//            }
        }

        r2.readLine();
        ArrayList<Integer> genMapPos = new ArrayList<Integer>();
        ArrayList<Double> genMapv1 = new ArrayList<Double>();
        while ((line = r2.readLine()) != null) {
            String[] tmp1 = line.split(" ");
            genMapPos.add(Integer.valueOf(tmp1[0]));
            genMapv1.add(Double.valueOf(tmp1[1]));
        }
        double mb = 100000000.0;
        int jstart = 0;
        for (int i = 1; i < pValues.size() - 1; i++) { 
            double recomrate = 0;

            if (pValues.get(i) >= genMapPos.get(genMapPos.size() - 1)) {
                recomrate = genMapv1.get(genMapv1.size() - 2) / mb;
            } else {
                for (int j = jstart; j < genMapPos.size() - 1; j++) {
                    //System.out.println( String.valueOf(pValues.get(i))+" "+String.valueOf(genMapPos.get(j))+" "+String.valueOf(genMapPos.get(j+1))      );
                    if (pValues.get(i) >= genMapPos.get(j) && pValues.get(i) < genMapPos.get(j + 1) && pValues.get(i + 1) <= genMapPos.get(j + 1)) {
                        recomrate = (double) genMapv1.get(j) / mb;

                        jstart = j;
                        break;
                    } else if (pValues.get(i) >= genMapPos.get(j) && pValues.get(i) < genMapPos.get(j + 1) && pValues.get(i + 1) > genMapPos.get(j + 1)) {
                        double recomCurrent = genMapv1.get(j) * ((double) genMapPos.get(j + 1) - (double) pValues.get(i));
                        int endspot = j + 1;

                        if (endspot == genMapPos.size() - 1) {
                            recomCurrent = recomCurrent + genMapv1.get(j) * ((double) pValues.get(i + 1) - (double) genMapPos.get(j + 1));
                            recomrate = (recomCurrent / ((double) pValues.get(i + 1) - (double) pValues.get(i))) / mb;
                            break;
                        }

                        while (pValues.get(i + 1) > genMapPos.get(endspot + 1)) {
                            recomCurrent = recomCurrent + genMapv1.get(endspot) * ((double) genMapPos.get(endspot + 1) - (double) genMapPos.get(endspot));
                            endspot = endspot + 1;
                            if (endspot == (genMapPos.size() - 1)) {
                                break;
                            }
                        }
                        recomCurrent = recomCurrent + genMapv1.get(endspot) * ((double) pValues.get(i + 1) - (double) genMapPos.get(endspot));
                        recomrate = (recomCurrent / ((double) pValues.get(i + 1) - pValues.get(i))) / mb;
                        jstart = j;
                    }

                }
            }
            if (recomrate == 0) {
                //System.out.println("Warning: zero recombination rate at SNP "+i);
                if ((pValues.get(i) < genMapPos.get(0)) && (pValues.get(i + 1) <= genMapPos.get(1))) {
                    recomrate = genMapv1.get(0) / mb;
                } else if ((pValues.get(i) < genMapPos.get(0)) && (pValues.get(i + 1) > genMapPos.get(1))) {
                    double recomCurrent = genMapv1.get(0) * ((double) genMapv1.get(1) - (double) pValues.get(i));
                    int endspot = 1;
                    while (pValues.get(i + 1) > genMapPos.get(endspot + 1)) {
                        recomCurrent = recomCurrent + genMapv1.get(endspot) * ((double) genMapPos.get(endspot + 1) - (double) genMapPos.get(endspot));
                        endspot = endspot + 1;
                        if (endspot == (genMapPos.size() - 1)) {
                            break;
                        }
                    }
                    recomCurrent = recomCurrent + genMapv1.get(endspot) * ((double) pValues.get(i + 1) - (double) genMapPos.get(endspot));
                    recomrate = (recomCurrent / ((double) pValues.get(i + 1) - (double) pValues.get(i))) / mb;
                }
            }
            wr2.write(String.valueOf(pValues.get(i)) + " " + String.format("%.15f", (recomrate)));
            wr2.newLine();
        }
        wr2.write(String.valueOf(pValues.get(pValues.size() - 1)) + " " + String.valueOf(0));
        wr2.newLine();
        wr2.close();

        int lin = (lineLength - 5);
        BufferedWriter wr1 = new BufferedWriter(new FileWriter(new File(outfile)));

        //wr1.write("0");
        //wr1.newLine();
        wr1.write(String.valueOf(lin));
        wr1.newLine();
        wr1.write(String.valueOf(count));
        wr1.newLine();

        wr1.write("P");
        for (int i = 1; i < pValues.size(); i++) {
            wr1.write(" " + String.valueOf(pValues.get(i)));
        }
        wr1.newLine();

        for (int i = 0; i < lin; i++) {
            for (BitSet set : sets) {
                if (set.get(i)) {
                    wr1.write("1");
                } else {
                    wr1.write("0");
                }
            }
            wr1.newLine();
        }
        wr1.close();
    }

    public static void withBitSetExperimenting(String filename, String outfile, String geneticMap, String recombOut) throws FileNotFoundException, IOException {
        BufferedReader r1 = new BufferedReader(new FileReader(new File(filename)));
        
        String line;
        //String[] tmp;
        int count = 0;
        //ArrayList<Integer> pValues = new ArrayList<Integer>();
        ArrayList<BitSet> sets = new ArrayList<BitSet>();
        int lineLength = 0;
        //pValues.add(-1);

        Runtime runtime = Runtime.getRuntime();

        NumberFormat format = NumberFormat.getInstance();

        //StringBuilder sb = new StringBuilder();
//        while ((line = r1.readLine()) != null) {
//            String[] tmp = line.split(" ");
//            lineLength = tmp.length;    
////            BitSet tk = new BitSet(90);
////            for (int i = 5; i < lineLength; i++) {
////                if (tmp[i].equals("1")) {
////                    tk.set(i - 5);
////                }
////            }
////            sets.add(tk);
//            count++;
//        }
        
        
        for (int i = 0; i < 1085830; i++) {
            BitSet tk = new BitSet(90);
            tk.set(0, 80);
            sets.add(tk);
        }

        //r1.close();

        System.out.println(count);
        int lin = (lineLength - 5);

        long maxMemory = runtime.maxMemory();
        long allocatedMemory = runtime.totalMemory();
        long freeMemory = runtime.freeMemory();

        System.out.println("free memory: " + format.format(freeMemory / 1024) + "<br/>");
        System.out.println("allocated memory: " + format.format(allocatedMemory / 1024) + "<br/>");
        System.out.println("max memory: " + format.format(maxMemory / 1024) + "<br/>");
        System.out.println("total free memory: " + format.format((freeMemory + (maxMemory - allocatedMemory))));

        for (int i = 0; i < lin; i++) {
            for (BitSet set : sets) {
                if (set.get(i)) {
                    //wr1.write("1");
                } else {
                    //wr1.write("0");
                }
            }

        }

    }
    
    public static void diagnostics(String file1, String file2) throws IOException
    {       
        BufferedReader r1 = new BufferedReader(new FileReader(new File(file1)));
        BufferedReader r2 = new BufferedReader(new FileReader(new File(file2)));
        String line1, line2;

        int count = 0;
        while ((line1 = r1.readLine()) != null) {
            line2 = r2.readLine();
            if (!line1.equals(line2)) {
                //System.out.println("mismatch");
                //System.out.println(line1.hashCode());
                //System.out.println(line2.hashCode());
                //System.out.println();
                //System.out.println("line :" + count);
                //System.out.println(line1 + " " + line2);

                double d1 = Double.valueOf(line1.split(" ")[1]);
                double d2 = Double.valueOf(line2.split(" ")[1]);

                if (Math.abs(d1 - d2) > 0.000000000000002) {
                    System.out.println("line :" + count);
                    System.out.println(line1 + " " + line2);
                    System.out.println(Math.abs(d1 - d2));
                }

                //System.out.println(line2.length());
//                String[] tmp1 = line1.split(" ");
//                String[] tmp2 = line2.split(" ");
//                System.out.println("elements: " + tmp1.length + " vs " + tmp2.length);
//                for (int i = 0; i < tmp1.length; i++) {
//                    if (!tmp1[i].equals(tmp2[i])) {
//                        System.out.println(tmp1[i] + " vs " + tmp2[i]);
//                    }
//                }
                //System.out.println("lengths "+tmp1.length+" "+tmp2.length);
//                for(int i=0;i<line1.length();i++)
//                {
//                    if(line1.charAt(i)!=line2.charAt(i))
//                    {
//                        System.out.println(line1.charAt(i)+" vs "+line2.charAt(i)+" @ "+i);
//                    }
//                }
            }
            count++;

        }
        byte[] f1 = Files.readAllBytes(new File("/home/thanos/chr21_3pobs.output").toPath());
        byte[] f2 = Files.readAllBytes(new File("/home/thanos/GOODchr21_3pobs.output").toPath());
        System.out.println(Arrays.equals(f1, f2));

        /**
         * * Diagnostics END **
         */
    }
    

    public static void main(String[] args) throws FileNotFoundException, IOException {

        //String flagZero = "1";
        //args[0]="/home/thanos/chr21_3pobs.haps";
        //args[1]="/home/thanos/chr21_3pobs.output";// args[1];
        //args[2]="/home/thanos/genetic_map_chr21_combined_b37.txt";
        
        String filename = args[0];
        String outfile = args[1];
        String geneticMap = args[2];
        String recomb = args[3];
        
        long t1 = System.currentTimeMillis();
        withBitSet(filename, outfile, geneticMap, recomb);
        System.out.println("Time : " + (System.currentTimeMillis() - t1));


    }

}
