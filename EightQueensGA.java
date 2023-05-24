import java.util.Random;
import java.util.Scanner;

/**
 * solves eight queens problem using genetic algorithms
 * 
 * @author Vance Spears
 * @version 3/29/23
 */
public class EightQueensGA {

   /** the size of the chess board */
   public static final int BOARD_SIZE = 8;

   /** the number of chromosomes in each generation */
   public static final int NUM_CHROMS = 20;

   /** the number of generations to produce */
   public static final int NUM_SIMS = 10000;

   /** the number of parents to reproduce children with */
   public static final int NUM_PARENTS = 3;

   /** the number of children to produce each generation */
   public static final int NUM_CHILDREN = 2;

   /** the number of gene mutations to perform on each generation's children */
   public static final int MUTATE_COUNT = 1;

   /** the percentage increment at which the simulation should pause */
   public static final int PAUSE_PERCENT = 10;

   /** an array of chromosomes representing the total population */
   public EightQueensChrom[] myChroms;

   /** a Random object for random number generation */
   public Random random = new Random();

   /** a Scanner object for reading user input */
   public Scanner reader = new Scanner(System.in);

   /**
    * creates a starting population of random chromosomes
    * 
    * @param num_chroms number of chromosomes in population
    */
   public void generate_initial_population(int num_chroms) {
      myChroms = new EightQueensChrom[num_chroms];
      for (int i = 0; i < num_chroms; i++) {
         int[] initialGenes = new int[BOARD_SIZE];
         for (int j = 0; j < BOARD_SIZE; j++) {
            initialGenes[j] = random.nextInt(BOARD_SIZE);
         }
         myChroms[i] = new EightQueensChrom(initialGenes);
      }
   }

   /**
    * creates the roulette wheel for selecting parents
    * 
    * @return array containing scaled selection probabilities based on fitness
    *         (total = 1.0)
    */
   public double[] set_probabilities_of_population() {
      double[] probabilities = new double[NUM_CHROMS];
      double[] reciprocalFitnesses = new double[NUM_CHROMS];
      double totalFitness = 0;

      for (int i = 0; i < NUM_CHROMS; i++) {
         double reciprocalFitness = 1 / myChroms[i].getFitness();
         reciprocalFitnesses[i] = reciprocalFitness;
         totalFitness += reciprocalFitness;
      }

      for (int i = 0; i < NUM_CHROMS; i++) {
         probabilities[i] = reciprocalFitnesses[i] / totalFitness;
      }

      return probabilities;
   }

   /**
    * picks random parents to use for creating offspring
    * 
    * @param number_of_selections number of parents to select
    * @return array containing selected parents
    */
   public EightQueensChrom[] roulette_wheel_selection(int number_of_selections) {
      EightQueensChrom[] parents = new EightQueensChrom[number_of_selections];
      double[] probabilities = set_probabilities_of_population();

      for (int i = 0; i < number_of_selections; i++) {
         double randomProbability = Math.random();
         boolean hasFoundParent = false;
         for (int j = 0; j < probabilities.length; j++) {
            if (!hasFoundParent) {
               if (randomProbability < probabilities[j]) {
                  parents[i] = myChroms[j];
                  hasFoundParent = true;
               } else {
                  randomProbability -= probabilities[j];
               }
            }
         }
         if (!hasFoundParent) {
            parents[i] = myChroms[0];
         }
      }

      return parents;
   }

   /**
    * produces children from the selected parents using one-point crossover
    * 
    * @param the_chosen   the array of selected parents
    * @param num_children the number of offspring to produce
    * @return the array of offspring
    */
   public EightQueensChrom[] reproduce_children(EightQueensChrom[] the_chosen, int num_children) {
      EightQueensChrom[] children = new EightQueensChrom[num_children];

      for (int i = 0; i < num_children; i++) {
         int firstRandom = 0;
         int secondRandom = 0;
         boolean shouldGenerateRandom = true;
         while (shouldGenerateRandom) {
            firstRandom = random.nextInt(the_chosen.length);
            secondRandom = random.nextInt(the_chosen.length);
            if (firstRandom != secondRandom) {
               shouldGenerateRandom = false;
            }
         }

         EightQueensChrom firstParent = the_chosen[firstRandom];
         EightQueensChrom secondParent = the_chosen[secondRandom];
         children[i] = crossWithTwoPoints(firstParent, secondParent);
      }

      return children;
   }

   /**
    * performs one-point crossover
    * 
    * @param firstParent  the first parent chromosome
    * @param secondParent the second parent chromosome
    * @return the child chromosome
    */
   EightQueensChrom crossWithOnePoint(EightQueensChrom firstParent, EightQueensChrom secondParent) {
      EightQueensChrom child = new EightQueensChrom();

      int cutoffIndex = random.nextInt(BOARD_SIZE);
      for (int i = 0; i < BOARD_SIZE; i++) {
         if (i <= cutoffIndex) {
            child.genes[i] = firstParent.genes[i];
         } else {
            child.genes[i] = secondParent.genes[i];
         }
      }

      return child;
   }

   /**
    * performs two-point crossover
    * 
    * @param firstParent  the first parent chromosome
    * @param secondParent the second parent chromosome
    * @return the child chromosome
    */
   EightQueensChrom crossWithTwoPoints(EightQueensChrom firstParent, EightQueensChrom secondParent) {
      int[] initGenes = new int[BOARD_SIZE];
      for (int i = 0; i < BOARD_SIZE; i++) {
         initGenes[i] = -1;
      }

      EightQueensChrom child = new EightQueensChrom(initGenes);

      int firstCutIndex = random.nextInt(BOARD_SIZE);
      int secondCutIndex = random.nextInt(BOARD_SIZE);
      int minCutIndex = Math.min(firstCutIndex, secondCutIndex);
      int maxCutIndex = Math.max(firstCutIndex, secondCutIndex);

      for (int i = minCutIndex; i <= maxCutIndex; i++) {
         child.genes[i] = secondParent.genes[i];
      }

      for (int i = 0; i < firstParent.genes.length; i++) {
         if (child.genes[i] == -1) {
            child.genes[i] = firstParent.genes[i];
         }
      }

      return child;
   }

   /**
    * mutates each offspring using single-point mutation (MUTATE_COUNT times)
    * 
    * @param the_children the array of children to mutate
    * @return the array of mutated children
    */
   public EightQueensChrom[] mutate_children(EightQueensChrom[] the_children) {
      for (int i = 0; i < the_children.length; i++) {
         for (int j = 0; j < MUTATE_COUNT; j++) {
            the_children[i].genes[random.nextInt(BOARD_SIZE)] = random.nextInt(BOARD_SIZE);
         }
      }

      return the_children;
   }

   /**
    * adds given children to the population by removing the weakest population
    * members
    * 
    * @param the_children array of new chromosomes to add to myChroms
    */
   public void merge_population_and_children(EightQueensChrom[] the_children) {
      sortChroms();
      for (int i = 0; i < the_children.length; i++) {
         myChroms[myChroms.length - i - 1] = the_children[i];
      }
   }

   /**
    * sorts myChroms in increasing order based on getFitness() values
    */
   public void sortChroms() {
      // using insertion sort
      for (int i = 1; i < myChroms.length; i++) {
         EightQueensChrom key = myChroms[i];
         int j = i;
         while (j > 0 && myChroms[j - 1].getFitness() > key.getFitness()) {
            myChroms[j] = myChroms[j - 1];
            j--;
         }

         myChroms[j] = key;
      }
   }

   /**
    * prints debug information and pauses program until enter is pressed
    * 
    * @param generation the current generation number
    */
   public void pause(int generation) {
      System.out.println();
      System.out.println("Population Statistics (gen " + generation
            + " of " + NUM_SIMS + "):");
      System.out.println("----------------------------------------");
      for (int i = 0; i < myChroms.length; i++) {
         System.out.println("chrom " + (i + 1) + ": " + myChroms[i]
               + ", fitness: " + myChroms[i].getFitness());
      }

      reader.nextLine();
   }

   /**
    * runs the ga simulation
    * 
    * @param population_size       number of chromosomes in population
    * @param number_of_generations number of generations to simulate
    * @return the fittest chromosome
    */
   public EightQueensChrom run_ga(int population_size, int number_of_generations) {
      double bestGlobalFitness = 0;
      generate_initial_population(population_size);

      for (int generation = 0; generation < number_of_generations; generation++) {
         if (100 * (double) generation / number_of_generations % PAUSE_PERCENT == 0) {
            pause(generation);
         }

         double currentBestFitness = 0;
         for (EightQueensChrom chrom : myChroms) {
            currentBestFitness += chrom.getFitness();
         }

         if (currentBestFitness > bestGlobalFitness) {
            bestGlobalFitness = currentBestFitness;
         }

         EightQueensChrom[] the_chosen = roulette_wheel_selection(NUM_PARENTS);
         EightQueensChrom[] the_children = reproduce_children(the_chosen, NUM_CHILDREN);
         the_children = mutate_children(the_children);
         merge_population_and_children(the_children);
      }

      sortChroms();
      return myChroms[0];
   }

   /**
    * runs the simulation
    */
   public void run() {
      EightQueensChrom best = run_ga(NUM_CHROMS, NUM_SIMS);
      System.out.println();
      System.out.println("Winner = " + best + " (fitness = " + best.getFitness()
            + ")");
      reader.close();
   }

   public static void main(String[] args) {
      EightQueensGA test = new EightQueensGA();
      test.run();
   }

   /**
    * a chromosome representation of the eight queens problem
    */
   private class EightQueensChrom {
      public int[] genes;

      /**
       * default constructor for EightQueensChrom
       */
      public EightQueensChrom() {
         genes = new int[BOARD_SIZE];
      }

      /**
       * an EightQueensChrom constructor for given gene values
       * 
       * @param geneVals the gene values
       */
      public EightQueensChrom(int[] geneVals) {
         this();
         for (int idx = 0; idx < genes.length; idx++) {
            genes[idx] = geneVals[idx];
         }
      }

      /**
       * gets fitness score of chromosome
       * 
       * @return number of queens each queen attacks
       */
      public double getFitness() {
         int fitnessScore = 0;

         for (int i = 0; i < genes.length; i++) {
            for (int j = 0; j < genes.length; j++) {
               if (i != j) {
                  // check if other queen is in row
                  if (genes[i] == genes[j]) {
                     fitnessScore++;
                  }
                  // check if other queen is in diagonal
                  if (genes[i] == genes[j] + (j - i)
                        || genes[i] == genes[j] - (j - i)) {
                     fitnessScore++;
                  }
               }
            }
         }

         return fitnessScore;
      }

      /**
       * gets string representation of chromosome
       * 
       * @return genes as a string
       */
      public String toString() {
         String out = "";

         for (int idx = 0; idx < genes.length; idx++) {
            out += genes[idx];
         }
         return out;
      }

   }
}