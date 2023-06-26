using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.IO;
using System.Numerics;
using System.Windows.Forms;


class Program
{
    public class GraphForm : Form
    {
        private List<Location> path;
        private double minX, maxX, minY, maxY;
        private const int Padding = 50;

        public GraphForm(List<Location> path)
        {
            this.path = path;

            // Find the minimum and maximum values of X and Y coordinates
            minX = maxX = path[0].X;
            minY = maxY = path[0].Y;

            foreach (Location location in path)
            {
                minX = Math.Min(minX, location.X);
                maxX = Math.Max(maxX, location.X);
                minY = Math.Min(minY, location.Y);
                maxY = Math.Max(maxY, location.Y);
            }

            // Subscribe to the Paint event to draw the graph
            Paint += GraphForm_Paint;
        }

        private void GraphForm_Paint(object sender, PaintEventArgs e)
        {
            // Create a Graphics object from the PaintEventArgs
            Graphics g = e.Graphics;

            // Calculate the scaling factors
            double xScale = (ClientSize.Width - 2 * Padding) / (maxX - minX);
            double yScale = (ClientSize.Height - 2 * Padding) / (maxY - minY);

            // Set up the pen for drawing the path
            using (Pen pen = new Pen(Color.Red, 2))
            {
                // Draw the path and calculate the total length
                double totalLength = 0;

                // Draw the path
                for (int i = 0; i < path.Count - 1; i++)
                {
                    Location start = path[i];
                    Location end = path[i + 1];

                    // Scale the coordinates
                    int startX = (int)((start.X - minX) * xScale) + Padding;
                    int startY = (int)((start.Y - minY) * yScale) + Padding;
                    int endX = (int)((end.X - minX) * xScale) + Padding;
                    int endY = (int)((end.Y - minY) * yScale) + Padding;

                    // Draw the line
                    g.DrawLine(pen, startX, startY, endX, endY);

                    // Calculate the midpoint for placing the length and name
                    int midX = (startX + endX) / 2;
                    int midY = (startY + endY) / 2;

                    // Calculate the length of the edge
                    double length = FindDistance(start, end);
                    totalLength += length;

                    // Draw the length
                    string lengthText = length.ToString("0.##");
                    g.DrawString(lengthText, Font, Brushes.Black, midX, midY);

                    // Draw the name of each point
                    g.DrawString(start.Name, Font, Brushes.Black, startX, startY);
                    g.DrawString(end.Name, Font, Brushes.Black, endX, endY);


                    


                }
                // Draw the total length at the bottom of the graph
                string totalLengthText = "Total Length: " + string.Format("{0:0.00}", totalLength);
                g.DrawString(totalLengthText, Font, Brushes.Black, 10, ClientSize.Height - 30);
            }
        }
        
       
        
    }



    public class Location
    {
        public string Name { get; set; }
        public string Id { get; set; }
        public double X { get; set; }
        public double Y { get; set; }

        public Location(string name, string id, double x, double y)
        {
            Name = name;
            Id = id;
            X = x;
            Y = y;
        }
    }
    private static double CalculateTotalCost(List<Location> path)
    {
        double totalCost = 0;

        for (int i = 0; i < path.Count - 1; i++)
        {
            Location currentLocation = path[i];
            Location nextLocation = path[i + 1];
            double distance = FindDistance(currentLocation, nextLocation);
            totalCost += distance;
        }

        return totalCost;
    }
    private static double FindDistance(Location l1, Location l2)
    {
        return Math.Sqrt(Math.Pow(l2.X - l1.X, 2) + Math.Pow(l2.Y - l1.Y, 2));
    }

    public static Location CreateLocation(string line)
    {
        var values = line.Split(',');
        var name = values[0];
        string id = (values[1]);
        var x = double.Parse(values[2]);
        var y = double.Parse(values[3]);

        return new Location(name, id, x, y);
    }

    public static List<Location> ReadLocationsFromFile(string path)
    {
        var locations = new List<Location>();

        using (var reader = new StreamReader(path))
        {

            // read each line and create a Location object
            while (!reader.EndOfStream)
            {
                var line = reader.ReadLine();
                var location = CreateLocation(line);
                locations.Add(location);
            }
        }

        return locations;
    }




    static void Main(string[] args)
    {
        var path = "data.csv";
        var locations = ReadLocationsFromFile(path);
        double[,] distances = CalculateDistances(locations);
        //List<Location> shortestPath = first.TSP_Implement(distances, 0, locations);
        int n = locations.Count;
        List<Location> shortestPath;
        List<Location> set1 = new List<Location>();
        List<Location> set2 = new List<Location>();
        List<Location> set3 = new List<Location>();
        List<Location> sp1 = null;
        List<Location> sp2 = null;
        List<Location> sp3 = null;
        for (int i = 0; i < n; i++)
        {
            if (i < n / 3) { set1.Add(locations[i]); }
            else if (i >= n / 3 && i < 2 * n / 3) { set2.Add(locations[i]); }
            else { set3.Add(locations[i]); }
        }
        // TimeSpan timeout = TimeSpan.FromSeconds(10);
        //CancellationTokenSource cancellationTokenSource = new CancellationTokenSource();
        //CancellationToken cancellationToken = cancellationTokenSource.Token;
        //Task.Run(() =>
        //{

        {
            Console.WriteLine("--Pirmas--");
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();

            sp1 = first.TSP_Implement(distances, 0, set1, locations);
            sp2 = first.TSP_Implement(distances, 0, set2, locations);
            sp3 = first.TSP_Implement(distances, 0, set3, locations);
            
            // });
            stopwatch.Stop();
            Application.Run(new GraphForm(sp1));
            Application.Run(new GraphForm(sp2));
            Application.Run(new GraphForm(sp3));
            //TimeSpan elapsed = stopwatch.Elapsed;
            //double secondsPassed = elapsed.TotalSeconds;
            //Console.WriteLine($"Seconds passed: {secondsPassed}");
            long elapsedMilliseconds = stopwatch.ElapsedMilliseconds;

            Console.WriteLine("Elapsed Time: {0} ms", elapsedMilliseconds);

            // Wait for the specified duration and then cancel the task
            //Task.Delay(timeout).ContinueWith(_ => cancellationTokenSource.Cancel());
            /*
            
            Console.WriteLine("     --1");
            if (sp1 != null)
            {
                for (int i = 0; i < sp1.Count(); i++)
                {
                    Location location = sp1[i];
                    if (i != 0)
                    {
                        { Console.WriteLine($"Dist: {FindDistance(sp1[i - 1], sp1[i])} ,Name: {location.Name},  X: {location.X}, Y: {location.Y}"); }
                    }
                    else
                    { Console.WriteLine($"Name: {location.Name},  X: {location.X}, Y: {location.Y}"); }
                }
            }
            else
            {
                Console.WriteLine("Could not find a valid path that visits all locations exactly once.");
            }
            Console.WriteLine("     --2");
            if (sp2 != null)
            {
                foreach (var location in sp2)
                {
                    Console.WriteLine($"Name: {location.Name}, ID: {location.Id}, X: {location.X}, Y: {location.Y}");
                }
            }
            else
            {
                Console.WriteLine("Could not find a valid path that visits all locations exactly once.");
            }
            Console.WriteLine("     --3");
            if (sp3 != null)
            {
                foreach (var location in sp3)
                {
                    Console.WriteLine($"Name: {location.Name}, ID: {location.Id}, X: {location.X}, Y: {location.Y}");
                }
            }
            else
            {
                Console.WriteLine("Could not find a valid path that visits all locations exactly once.");
            } */
        }
        // Console.WriteLine(CalculateTotalCost(shortestPath));
        {
            Console.WriteLine();
            Console.WriteLine("--Antras--");
            /*shortestPath = second.GreedyTSP(distances, 0, locations, locations);
            if (shortestPath != null)
            {
                foreach (var location in shortestPath)
                {
                    Console.WriteLine($"Name: {location.Name}, ID: {location.Id}, X: {location.X}, Y: {location.Y}");
                }
            }
            else
            {
                Console.WriteLine("Could not find a valid path that visits all locations exactly once.");
            } */
            // Console.WriteLine(CalculateTotalCost(shortestPath));
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();
            sp1 = second.GreedyTSP(distances, 0, set1, locations);
            sp2 = second.GreedyTSP(distances, 0, set2, locations);
            sp3 = second.GreedyTSP(distances, 0, set3, locations);

            stopwatch.Stop();
            Application.Run(new GraphForm(sp1));
            Application.Run(new GraphForm(sp2));
            Application.Run(new GraphForm(sp3));
            long elapsedMilliseconds = stopwatch.Elapsed.Milliseconds;
            Console.WriteLine("Elapsed Time: {0} ms", elapsedMilliseconds);
            /*
            

            Console.WriteLine("     --1");
            if (sp1 != null)
            {
                foreach (var location in sp1)
                {
                    Console.WriteLine($"Name: {location.Name}, ID: {location.Id}, X: {location.X}, Y: {location.Y}");
                }
            }
            else
            {
                Console.WriteLine("Could not find a valid path that visits all locations exactly once.");
            }
            Console.WriteLine("     --2");
            if (sp2 != null)
            {
                foreach (var location in sp2)
                {
                    Console.WriteLine($"Name: {location.Name}, ID: {location.Id}, X: {location.X}, Y: {location.Y}");
                }
            }
            else
            {
                Console.WriteLine("Could not find a valid path that visits all locations exactly once.");
            }
            Console.WriteLine("     --3");
            if (sp3 != null)
            {
                foreach (var location in sp3)
                {
                    Console.WriteLine($"Name: {location.Name}, ID: {location.Id}, X: {location.X}, Y: {location.Y}");
                }
            }
            else
            {
                Console.WriteLine("Could not find a valid path that visits all locations exactly once.");
            }
            */
        }
        { 
            Console.WriteLine();
            Console.WriteLine("--Trecias--");
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();
            sp1 = third.GeneticTSP(distances, 0, set1, 100, 100);
            sp2 = third.GeneticTSP(distances, 0, set2, 100, 100);
            sp3 = third.GeneticTSP(distances, 0, set3, 100, 100);
            //shortestPath = third.GeneticTSP(distances, 0, locations, 100, 100);
            stopwatch.Stop();
           // Application.Run(new GraphForm(sp1));
            //Application.Run(new GraphForm(sp2));
            //Application.Run(new GraphForm(sp3));
            long elapsedMilliseconds = stopwatch.ElapsedMilliseconds;
            Console.WriteLine("Elapsed Time: {0} ms", elapsedMilliseconds);
            /*if (shortestPath != null)
            {
                foreach (var location in shortestPath)
                {
                    Console.WriteLine($"Name: {location.Name}, ID: {location.Id}, X: {location.X}, Y: {location.Y}");
                }
            }
            else
            {
                Console.WriteLine("Could not find a valid path that visits all locations exactly once.");
            }*/
        }
    }
    static double[,] CalculateDistances(List<Location> places)
    {
        var distances = new double[places.Count, places.Count];

        for (int i = 0; i < places.Count; i++)
        {
            for (int j = i + 1; j < places.Count; j++)
            {
                var distance = FindDistance(places[i], places[j]);
                distances[i, j] = distance;
                distances[j, i] = distance;
            }
        }

        return distances;
    }
    public class first
    {
        //keliaujancio pirklio problem sprend.
        public static List<Location> TSP_Implement(double[,] adjMatrix, int start, List<Location> locations, List<Location> original)
        {
            //locations.Insert(0, original[start]);
            int V = locations.Count;
            var cities = new List<int>();
            for (int i = 0; i < V; i++)
            {
                if (i != start)
                {
                    cities.Add(i);
                }
            }

            double minDistance = double.MaxValue;
            List<Location> shortestPath = null;

            do
            {
                double currDistance = 0;
                int k = start;

                var path = new List<Location> { locations[start] };

                for (int i = 0; i < cities.Count; i++)
                {
                    currDistance += adjMatrix[k, cities[i]];
                    k = cities[i];
                    path.Add(locations[cities[i]]);
                }

                currDistance += adjMatrix[k, start];
                path.Add(locations[start]); // Add the starting location to the end of the path

                if (currDistance < minDistance)
                {
                    minDistance = currDistance;
                    shortestPath = path;
                }
            } while (NextPermutation(cities));

            return shortestPath;
        }
        //permutacijos, jos sugeneruota visus galimus maršrutus ir paskaičiuotų kiekvieno maršruto bendrą atstumą. 
        //Pradedant nuo pradinio miesto, jis eina per visus likusius miestus ir skaičiuoja atstumą tarp kiekvieno miesto ir jo sekančio miesto maršrute.
        //Metodas vykdo permutacijų operaciją, kol visi galimi maršrutai yra išbandyti.
        private static bool NextPermutation(List<int> cities)
        {
            int i = cities.Count - 2;
            while (i >= 0 && cities[i] >= cities[i + 1])
            {
                i--;
            }

            if (i < 0)
            {
                return false;
            }

            int j = cities.Count - 1;
            while (cities[j] <= cities[i])
            {
                j--;
            }

            Swap(cities, i, j);
            Reverse(cities, i + 1, cities.Count - 1);
            return true;
        }
        private static void Swap(List<int> cities, int i, int j)
        {
            int temp = cities[i];
            cities[i] = cities[j];
            cities[j] = temp;
        }

        private static void Reverse(List<int> cities, int start, int end)
        {
            while (start < end)
            {
                Swap(cities, start, end);
                start++;
                end--;
            }
        }


    }
    public class second
    {
        public static List<Location> GreedyTSP(double[,] adjMatrix, int start, List<Location> locations, List<Location> original)
        {
            int n = locations.Count;
            bool[] visited = new bool[n];
            List<Location> path = new List<Location>();
            //locations.Insert(0, original[start]);
            path.Add(locations[start]);
            visited[start] = true;

            while (path.Count < n)
            {
                int lastCityIndex = path.Count - 1;
                double minDistance = double.MaxValue;
                int closestCityIndex = -1;

                for (int i = 0; i < n; i++)
                {
                    if (!visited[i] && adjMatrix[lastCityIndex, i] < minDistance)
                    {
                        minDistance = adjMatrix[lastCityIndex, i];
                        closestCityIndex = i;
                    }
                }

                if (closestCityIndex != -1)
                {
                    path.Add(locations[closestCityIndex]);
                    visited[closestCityIndex] = true;
                }
            }

            // path.Add(original[start]); // Add the starting location to complete the cycle
            path.Add(locations[start]);
            return path;
        }
    }
    public class third
    {
        static Random random = new Random();

        public static List<Location> GeneticTSP(double[,] adjMatrix, int start, List<Location> locations, int populationSize, int maxGenerations)
        {
            int n = locations.Count;
            List<int[]> population = InitializePopulation(n, populationSize);
            List<Location> bestPath = null;
            double bestFitness = double.MaxValue;

            for (int generation = 0; generation < maxGenerations; generation++)
            {
                List<double> fitnessValues = CalculateFitness(adjMatrix, population, locations);
                int eliteIndex = GetEliteIndex(fitnessValues);
                double eliteFitness = fitnessValues[eliteIndex];

                if (eliteFitness < bestFitness)
                {
                    bestFitness = eliteFitness;
                    bestPath = ConvertToPath(population[eliteIndex], locations);
                }

                List<int[]> newPopulation = new List<int[]>();

                while (newPopulation.Count < populationSize)
                {
                    int[] parent1 = SelectParent(population, fitnessValues);
                    int[] parent2 = SelectParent(population, fitnessValues);

                    int[] child = Crossover(parent1, parent2);
                    child = Mutate(child);
                    newPopulation.Add(child);
                }

                population = newPopulation;
            }

            bestPath.Add(locations[start]); // Add the starting location to complete the cycle

            return bestPath;
        }

        static List<int[]> InitializePopulation(int n, int populationSize)
        {
            List<int[]> population = new List<int[]>();

            for (int i = 0; i < populationSize; i++)
            {
                int[] individual = Enumerable.Range(0, n).ToArray();
                ShuffleArray(individual);
                population.Add(individual);
            }

            return population;
        }

        static void ShuffleArray(int[] array)
        {
            int n = array.Length;

            for (int i = 0; i < n - 1; i++)
            {
                int j = random.Next(i, n);
                int temp = array[i];
                array[i] = array[j];
                array[j] = temp;
            }
        }

        static List<double> CalculateFitness(double[,] adjMatrix, List<int[]> population, List<Location> locations)
        {
            List<double> fitnessValues = new List<double>();

            foreach (int[] individual in population)
            {
                double distance = 0;

                for (int i = 0; i < individual.Length - 1; i++)
                {
                    int city1 = individual[i];
                    int city2 = individual[i + 1];
                    distance += adjMatrix[city1, city2];
                }

                distance += adjMatrix[individual[individual.Length - 1], individual[0]];
                fitnessValues.Add(1 / distance);
            }

            return fitnessValues;
        }

        static int GetEliteIndex(List<double> fitnessValues)
        {
            int eliteIndex = 0;
            double maxFitness = double.MinValue;

            for (int i = 0; i < fitnessValues.Count; i++)
            {
                if (fitnessValues[i] > maxFitness)
                {
                    maxFitness = fitnessValues[i];
                    eliteIndex = i;
                }
            }

            return eliteIndex;
        }

        static List<Location> ConvertToPath(int[] individual, List<Location> locations)
        {
            List<Location> path = new List<Location>();
            foreach (int cityIndex in individual)
            {
                path.Add(locations[cityIndex]);
            }

            return path;
        }


        static int[] SelectParent(List<int[]> population, List<double> fitnessValues)
        {
            int index = RouletteWheelSelection(fitnessValues);
            return population[index];
        }

        static int RouletteWheelSelection(List<double> fitnessValues)
        {
            double totalFitness = fitnessValues.Sum();
            double randomValue = random.NextDouble() * totalFitness;
            double partialSum = 0;
            for (int i = 0; i < fitnessValues.Count; i++)
            {
                partialSum += fitnessValues[i];

                if (partialSum >= randomValue)
                {
                    return i;
                }
            }

            return fitnessValues.Count - 1;
        }

        static int[] Crossover(int[] parent1, int[] parent2)
        {
            int n = parent1.Length;
            int[] child = new int[n];
            int startPos = random.Next(n);
            int endPos = random.Next(n);
            for (int i = 0; i < n; i++)
            {
                if (startPos < endPos && i > startPos && i < endPos)
                {
                    child[i] = parent1[i];
                }
                else if (startPos > endPos && !(i < startPos && i > endPos))
                {
                    child[i] = parent1[i];
                }
            }

            for (int i = 0; i < n; i++)
            {
                if (!child.Contains(parent2[i]))
                {
                    for (int j = 0; j < n; j++)
                    {
                        if (child[j] == 0)
                        {
                            child[j] = parent2[i];
                            break;
                        }
                    }
                }
            }

            return child;
        }
        static int[] Mutate(int[] individual)
        {
            int n = individual.Length;
            int mutationPoint1 = random.Next(n);
            int mutationPoint2 = random.Next(n);
            int temp = individual[mutationPoint1];
            individual[mutationPoint1] = individual[mutationPoint2];
            individual[mutationPoint2] = temp;
            return individual;
        }
    }
}