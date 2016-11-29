using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Diagnostics;


namespace TSP
{

    class ProblemAndSolver
    {

        private class TSPSolution
        {
            /// <summary>
            /// we use the representation [cityB,cityA,cityC] 
            /// to mean that cityB is the first city in the solution, cityA is the second, cityC is the third 
            /// and the edge from cityC to cityB is the final edge in the path.  
            /// You are, of course, free to use a different representation if it would be more convenient or efficient 
            /// for your data structure(s) and search algorithm. 
            /// </summary>
            public ArrayList
                Route;

            /// <summary>
            /// constructor
            /// </summary>
            /// <param name="iroute">a (hopefully) valid tour</param>
            public TSPSolution(ArrayList iroute)
            {
                Route = new ArrayList(iroute);
            }

            /// <summary>
            /// Compute the cost of the current route.  
            /// Note: This does not check that the route is complete.
            /// It assumes that the route passes from the last city back to the first city. 
            /// </summary>
            /// <returns></returns>
            public double costOfRoute()
            {
                // go through each edge in the route and add up the cost. 
                int x;
                City here;
                double cost = 0D;

                for (x = 0; x < Route.Count - 1; x++)
                {
                    here = Route[x] as City;
                    cost += here.costToGetTo(Route[x + 1] as City);
                }

                // go from the last city to the first. 
                here = Route[Route.Count - 1] as City;
                cost += here.costToGetTo(Route[0] as City);
                return cost;
            }
        }

        #region Private members 

        /// <summary>
        /// Default number of cities (unused -- to set defaults, change the values in the GUI form)
        /// </summary>
        // (This is no longer used -- to set default values, edit the form directly.  Open Form1.cs,
        // click on the Problem Size text box, go to the Properties window (lower right corner), 
        // and change the "Text" value.)
        private const int DEFAULT_SIZE = 25;

        /// <summary>
        /// Default time limit (unused -- to set defaults, change the values in the GUI form)
        /// </summary>
        // (This is no longer used -- to set default values, edit the form directly.  Open Form1.cs,
        // click on the Time text box, go to the Properties window (lower right corner), 
        // and change the "Text" value.)
        private const int TIME_LIMIT = 60;        //in seconds

        private const int CITY_ICON_SIZE = 5;


        // For normal and hard modes:
        // hard mode only
        private const double FRACTION_OF_PATHS_TO_REMOVE = 0.20;

        /// <summary>
        /// the cities in the current problem.
        /// </summary>
        private City[] Cities;
        /// <summary>
        /// a route through the current problem, useful as a temporary variable. 
        /// </summary>
        private ArrayList Route;
        /// <summary>
        /// best solution so far. 
        /// </summary>
        private TSPSolution bssf;

        /// <summary>
        /// how to color various things. 
        /// </summary>
        private Brush cityBrushStartStyle;
        private Brush cityBrushStyle;
        private Pen routePenStyle;


        /// <summary>
        /// keep track of the seed value so that the same sequence of problems can be 
        /// regenerated next time the generator is run. 
        /// </summary>
        private int _seed;
        /// <summary>
        /// number of cities to include in a problem. 
        /// </summary>
        private int _size;

        /// <summary>
        /// Difficulty level
        /// </summary>
        private HardMode.Modes _mode;

        /// <summary>
        /// random number generator. 
        /// </summary>
        private Random rnd;

        /// <summary>
        /// time limit in milliseconds for state space search
        /// can be used by any solver method to truncate the search and return the BSSF
        /// </summary>
        private int time_limit;
        #endregion

        #region Public members

        /// <summary>
        /// These three constants are used for convenience/clarity in populating and accessing the results array that is passed back to the calling Form
        /// </summary>
        public const int COST = 0;
        public const int TIME = 1;
        public const int COUNT = 2;

        public int Size
        {
            get { return _size; }
        }

        public int Seed
        {
            get { return _seed; }
        }
        #endregion

        #region Constructors
        public ProblemAndSolver()
        {
            this._seed = 1;
            rnd = new Random(1);
            this._size = DEFAULT_SIZE;
            this.time_limit = TIME_LIMIT * 1000;                  // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }

        public ProblemAndSolver(int seed)
        {
            this._seed = seed;
            rnd = new Random(seed);
            this._size = DEFAULT_SIZE;
            this.time_limit = TIME_LIMIT * 1000;                  // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }

        public ProblemAndSolver(int seed, int size)
        {
            this._seed = seed;
            this._size = size;
            rnd = new Random(seed);
            this.time_limit = TIME_LIMIT * 1000;                        // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }
        public ProblemAndSolver(int seed, int size, int time)
        {
            this._seed = seed;
            this._size = size;
            rnd = new Random(seed);
            this.time_limit = time * 1000;                        // time is entered in the GUI in seconds, but timer wants it in milliseconds

            this.resetData();
        }
        #endregion

        #region Private Methods

        /// <summary>
        /// Reset the problem instance.
        /// </summary>
        private void resetData()
        {

            Cities = new City[_size];
            Route = new ArrayList(_size);
            bssf = null;

            if (_mode == HardMode.Modes.Easy)
            {
                for (int i = 0; i < _size; i++)
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble());
            }
            else // Medium and hard
            {
                for (int i = 0; i < _size; i++)
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble(), rnd.NextDouble() * City.MAX_ELEVATION);
            }

            HardMode mm = new HardMode(this._mode, this.rnd, Cities);
            if (_mode == HardMode.Modes.Hard)
            {
                int edgesToRemove = (int)(_size * FRACTION_OF_PATHS_TO_REMOVE);
                mm.removePaths(edgesToRemove);
            }
            City.setModeManager(mm);

            cityBrushStyle = new SolidBrush(Color.Black);
            cityBrushStartStyle = new SolidBrush(Color.Red);
            routePenStyle = new Pen(Color.Blue, 1);
            routePenStyle.DashStyle = System.Drawing.Drawing2D.DashStyle.Solid;
        }

        #endregion

        #region Public Methods

        /// <summary>
        /// make a new problem with the given size.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size, HardMode.Modes mode)
        {
            this._size = size;
            this._mode = mode;
            resetData();
        }

        /// <summary>
        /// make a new problem with the given size, now including timelimit paremeter that was added to form.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size, HardMode.Modes mode, int timelimit)
        {
            this._size = size;
            this._mode = mode;
            this.time_limit = timelimit * 1000;                                   //convert seconds to milliseconds
            resetData();
        }

        /// <summary>
        /// return a copy of the cities in this problem. 
        /// </summary>
        /// <returns>array of cities</returns>
        public City[] GetCities()
        {
            City[] retCities = new City[Cities.Length];
            Array.Copy(Cities, retCities, Cities.Length);
            return retCities;
        }

        /// <summary>
        /// draw the cities in the problem.  if the bssf member is defined, then
        /// draw that too. 
        /// </summary>
        /// <param name="g">where to draw the stuff</param>
        public void Draw(Graphics g)
        {
            float width = g.VisibleClipBounds.Width - 45F;
            float height = g.VisibleClipBounds.Height - 45F;
            Font labelFont = new Font("Arial", 10);

            // Draw lines
            if (bssf != null)
            {
                // make a list of points. 
                Point[] ps = new Point[bssf.Route.Count];
                int index = 0;
                foreach (City c in bssf.Route)
                {
                    if (index < bssf.Route.Count - 1)
                        g.DrawString(" " + index + "(" + c.costToGetTo(bssf.Route[index + 1] as City) + ")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    else
                        g.DrawString(" " + index + "(" + c.costToGetTo(bssf.Route[0] as City) + ")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    ps[index++] = new Point((int)(c.X * width) + CITY_ICON_SIZE / 2, (int)(c.Y * height) + CITY_ICON_SIZE / 2);
                }

                if (ps.Length > 0)
                {
                    g.DrawLines(routePenStyle, ps);
                    g.FillEllipse(cityBrushStartStyle, (float)Cities[0].X * width - 1, (float)Cities[0].Y * height - 1, CITY_ICON_SIZE + 2, CITY_ICON_SIZE + 2);
                }

                // draw the last line. 
                g.DrawLine(routePenStyle, ps[0], ps[ps.Length - 1]);
            }

            // Draw city dots
            foreach (City c in Cities)
            {
                g.FillEllipse(cityBrushStyle, (float)c.X * width, (float)c.Y * height, CITY_ICON_SIZE, CITY_ICON_SIZE);
            }

        }

        /// <summary>
        ///  return the cost of the best solution so far. 
        /// </summary>
        /// <returns></returns>
        public double costOfBssf()
        {
            if (bssf != null)
                return (bssf.costOfRoute());
            else
                return -1D;
        }

        /// <summary>
        /// This is the entry point for the default solver
        /// which just finds a valid random tour 
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] defaultSolveProblem()
        {
            int i, swap, temp, count = 0;
            string[] results = new string[3];
            int[] perm = new int[Cities.Length];
            Route = new ArrayList();
            Random rnd = new Random();
            Stopwatch timer = new Stopwatch();

            timer.Start();

            do
            {
                for (i = 0; i < perm.Length; i++)                                 // create a random permutation template
                    perm[i] = i;
                for (i = 0; i < perm.Length; i++)
                {
                    swap = i;
                    while (swap == i)
                        swap = rnd.Next(0, Cities.Length);
                    temp = perm[i];
                    perm[i] = perm[swap];
                    perm[swap] = temp;
                }
                Route.Clear();
                for (i = 0; i < Cities.Length; i++)                            // Now build the route using the random permutation 
                {
                    Route.Add(Cities[perm[i]]);
                }
                bssf = new TSPSolution(Route);
                count++;
            } while (costOfBssf() == double.PositiveInfinity);                // until a valid route is found
            timer.Stop();

            results[COST] = costOfBssf().ToString();                          // load results array
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = count.ToString();

            return results;
        }

        //the id is just a unique identifier for each state_data to make it work with my implementation of the priority queue more easily
        int id = 0;
        /// <summary>
        /// performs a Branch and Bound search of the state space of partial tours
        /// stops when time limit expires and uses BSSF as solution
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] bBSolveProblem()
        {
            string[] results = new string[3];

            results[COST] = "not implemented";    // load results into array here, replacing these dummy values
            results[TIME] = "-1";
            results[COUNT] = "-1";

            int count = 0;

            Stopwatch timer = new Stopwatch();
            timer.Start();

            PriorityQueue queue = new PriorityQueue();
            queue.make_queue(Cities.Length * Cities.Length * Cities.Length * Cities.Length);

            //Here I run the greedy algorithm and use that BSSF result as my starting BSSF
            greedySolveProblem();
            Console.Out.WriteLine("Greedy BSSF: " + costOfBssf());

            matrix_and_bound initial_matrix = construct_initial_matrix();

            /*
            for (int j = 0; j < Cities.Length;j++)
            {
                for (int h = 0; h < Cities.Length; h++)
                {
                    Console.Write(initial_matrix.Matrix[j,h] + ",");
                }
                Console.WriteLine();
            }*/
            Console.WriteLine(initial_matrix.Lower_bound);
            //Console.WriteLine();


            //This loop will do my initial population of the queue by looping through each city and getting the lower bound
            for (int i = 0; i < Cities.Length; i++)
            {
                for (int k = 0; k < Cities.Length; k++)
                {
                    if (i != k)
                    {
                        //I need to get the lower bound of each reduced matrix and the bound
                        matrix_and_bound current = matrix_reduction(initial_matrix, i, k);

                        /*
                        for (int j = 0; j < Cities.Length; j++)
                        {
                            for (int h = 0; h < Cities.Length; h++)
                            {
                                Console.Write(current.Matrix[j, h] + ",");
                            }
                            Console.WriteLine();
                        }*/
                        //Console.WriteLine(current.Lower_bound);
                        //Console.WriteLine();

                        //If the lower bound is less than current bssf add to queue for later checking, otherwise ignore
                        if (current.Lower_bound < costOfBssf())
                        {

                            //need to create new state_data object with state data set
                            state_data data = new state_data();
                            //I guess depth doesn't matter to be exact so long as it's relative, so I'll keep this first one as 0
                            data.Depth = 0;
                            data.Mb = current;

                            data.add_city(Cities[i]);
                            data.add_city(Cities[k]);

                            data.add_city_index(i);
                            data.add_city_index(k);

                            data.set_priority();
                            //Console.Out.WriteLine("Set Priority " + data.Priority);
                            //I'm not sure this id is necessary but I'll have to see
                            data.Id = id;

                            queue.insert(data, id);
                            id++;
                        }
                    }
                }
            }


            //now run while queue is not empty and timer is less than 60
            while (timer.Elapsed.TotalSeconds < 60 && queue.Length > 0)
            {
                //pop off of queue and repeat above matrix reduction process
                state_data current = queue.delete_min();
                //if it's greater or equal to the current bssf I can just ignore it, this is the culling
                if (current.Mb.Lower_bound < costOfBssf())
                {
                    //Priority queue looks like it's working
                    //Console.Out.WriteLine(current.Mb.Lower_bound);

                    /*
                     * Now I need to matrix reduce each child of the current node, 
                     * see if it's still smaller than the current bssf if it's a leaf node
                     * I put it as the current bssf, otherwise I push it on the queue
                     */

                    for (int k = 0; k < Cities.Length; k++)
                    {
                        if (!current.City_list.Contains(k))
                        {
                            matrix_and_bound child = matrix_reduction(current.Mb, (int)current.City_list[current.City_list.Count - 1], k);

                            //Console.Out.WriteLine("here");

                            if (child.Lower_bound < costOfBssf())
                            {

                                //need to create new state_data object with state data set
                                state_data data = new state_data();
                                //I guess depth doesn't matter to be exact so long as it's relative, so I'll keep this first one as 0
                                data.Depth = current.Depth + 1;
                                data.Mb = child;

                                //The last value in the path is the current city
                                data.Path = current.Path;
                                data.add_city(Cities[k]);

                                data.City_list = current.City_list;
                                data.add_city_index(k);

                                data.set_priority();
                                //Console.Out.WriteLine("Set Priority " + data.Priority);
                                //I'm not sure this id is necessary but I'll have to see
                                data.Id = id;
                                id++;

                               // Console.Out.WriteLine("Intermediate Lower Bound " + data.Mb.Lower_bound);

                                if (data.City_list.Count < Cities.Length)
                                {
                                    queue.insert(data, id);
                                    
                                }
                                else if (data.Mb.Lower_bound < costOfBssf())//it's a leaf node and it's less than the current BSSF
                                {
                                    Console.Out.WriteLine("Current Lower Bound " + costOfBssf());
                                    Console.Out.WriteLine("Final Lower Bound " + data.Mb.Lower_bound);
                                    /*for (int j = 0; j < data.City_list.Count; j++)
                                    {
                                        Console.Out.WriteLine(data.City_list[j]);
                                    }*/

                                    
                                    bssf = new TSPSolution(data.Path);
                                    count++;
                                }
                                
                            }

                        }
                    }



                }
            }

            results[COST] = costOfBssf().ToString();
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = count.ToString();

            return results;
        }

        private matrix_and_bound construct_initial_matrix()
        {
            double[,] ret_matrix = new double[Cities.Length, Cities.Length];
            for (int i = 0; i < Cities.Length; i++)
            {
                for (int k = 0; k < Cities.Length; k++)
                {
                    if (i == k)
                    {
                        ret_matrix[i,k] = double.PositiveInfinity;
                    }
                    else
                    {
                        ret_matrix[i,k] = Cities[i].costToGetTo(Cities[k]);
                    }
                    
                }
            }

            double lower_bound = 0;

            /*
             *  Goes through all the rows and reduces them 
             */
            for (int i = 0; i < Cities.Length; i++)
            {
                double lowest = -1;
                for (int k = 0; k < Cities.Length; k++)
                {
                    if (!double.IsPositiveInfinity(ret_matrix[i,k]))
                    {
                        if (lowest == -1)
                        {
                            lowest = ret_matrix[i, k];
                        }
                        else if (lowest > ret_matrix[i, k])
                        {
                            lowest = ret_matrix[i, k];
                        }
                    }

                }
                lower_bound += lowest;
                //now subtract each value by the lowest
                for (int k = 0; k < Cities.Length; k++)
                {
                    if (!double.IsPositiveInfinity(ret_matrix[i, k]))
                    {
                        ret_matrix[i, k] = ret_matrix[i, k] - lowest;
                    }

                }
            }
            
            /*
             *  Goes through all the columns and reduces them 
             */
            for (int i = 0; i < Cities.Length; i++)
            {
                double lowest = -1;
                for (int k = 0; k < Cities.Length; k++)
                {
                    if (!double.IsPositiveInfinity(ret_matrix[k, i]))
                    {
                        if (lowest == -1)
                        {
                            lowest = ret_matrix[k, i];
                        }
                        else if (lowest > ret_matrix[k, i])
                        {
                            lowest = ret_matrix[k, i];
                        }
                    }

                }
                lower_bound += lowest;
                //now subtract each value by the lowest
                for (int k = 0; k < Cities.Length; k++)
                {
                    if (!double.IsPositiveInfinity(ret_matrix[k, i]))
                    {
                        ret_matrix[k, i] = ret_matrix[k, i] - lowest;
                    }

                }
            }


            return new matrix_and_bound(ret_matrix, lower_bound);
        }

        private matrix_and_bound matrix_reduction(matrix_and_bound matrix_in, int from, int to)
        {
            //create reduced matrix from the current city, will possibly need the array of available cities as well
            double[,] matrix = (double[,])matrix_in.Matrix.Clone();
            //Array.Copy(matrix_in.Matrix, matrix, Cities.Length); 
            double lower_bound = matrix_in.Lower_bound;

            if (! double.IsPositiveInfinity(matrix[from, to]))
            {
                lower_bound += matrix[from, to];
                matrix[from, to] = double.PositiveInfinity;

                //set row i and column k to infinity
                for (int h = 0; h < Cities.Length; h++)
                {
                    matrix[from, h] = double.PositiveInfinity;
                    matrix[h, to] = double.PositiveInfinity;
                }

                //now need to reduce the matrix again
                /*
                 *  Goes through all the rows and reduces them 
                 */
                for (int i = 0; i < Cities.Length; i++)
                {
                    double lowest = -1;
                    for (int k = 0; k < Cities.Length; k++)
                    {
                        if (!double.IsPositiveInfinity(matrix[i, k]))
                        {
                            if (lowest == -1)
                            {
                                lowest = matrix[i, k];
                            }
                            else if (lowest > matrix[i, k])
                            {
                                lowest = matrix[i, k];
                            }
                        }

                    }
                    if (lowest != -1)
                    {
                        lower_bound += lowest;
                        //now subtract each value by the lowest
                        for (int k = 0; k < Cities.Length; k++)
                        {
                            if (!double.IsPositiveInfinity(matrix[i, k]))
                            {
                                matrix[i, k] = matrix[i, k] - lowest;
                            }

                        }
                    }
                }

                /*
                 *  Goes through all the columns and reduces them 
                 */
                for (int i = 0; i < Cities.Length; i++)
                {
                    double lowest = -1;
                    for (int k = 0; k < Cities.Length; k++)
                    {
                        if (!double.IsPositiveInfinity(matrix[k, i]))
                        {
                            if (lowest == -1)
                            {
                                lowest = matrix[k, i];
                            }
                            else if (lowest > matrix[k, i])
                            {
                                lowest = matrix[k, i];
                            }
                        }

                    }
                    if (lowest != -1)
                    {
                        lower_bound += lowest;
                        //now subtract each value by the lowest
                        for (int k = 0; k < Cities.Length; k++)
                        {
                            if (!double.IsPositiveInfinity(matrix[k, i]))
                            {
                                matrix[k, i] = matrix[k, i] - lowest;
                            }

                        }
                    }
                }

            }
            else
            {
                lower_bound = double.PositiveInfinity;
            }



            return new matrix_and_bound(matrix, lower_bound);
        }

        /////////////////////////////////////////////////////////////////////////////////////////////
        // These additional solver methods will be implemented as part of the group project.
        ////////////////////////////////////////////////////////////////////////////////////////////

        /// <summary>
        /// finds the greedy tour starting from each city and keeps the best (valid) one
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] greedySolveProblem()
        {
            string[] results = new string[3];

            // TODO: Add your implementation for a greedy solver here.

            results[COST] = "not implemented";    // load results into array here, replacing these dummy values
            results[TIME] = "-1";
            results[COUNT] = "-1";
            int count = 0;

            TSPSolution solution = null;

            Stopwatch timer = new Stopwatch();

            timer.Start();
            //Outer loop for looking through each city
            for (int i = 0; i < Cities.Length; i++)
            {

                Route = new ArrayList();

                //This array will contain cities not visited already
                City[] current_cities = GetCities();

                //removes the city that is currently being looked at
                var temp = new List<City>(current_cities);
                temp.RemoveAt(i);
                City[] unchecked_cities = temp.ToArray();

                City cur_city = Cities[i];
                Route.Add(cur_city);

                //holds the shortest distance and the index of that city
                double shortest_distance = -1;
                int index = 0;

                do
                {
                    //Loops through the cities not visited yet to find smallest
                    for (int k = 0; k < unchecked_cities.Length; k++)
                    {
                        if (shortest_distance == -1 && !double.IsPositiveInfinity(cur_city.costToGetTo(unchecked_cities[k])))
                        {
                            shortest_distance = cur_city.costToGetTo(unchecked_cities[k]);
                        }
                        else if (shortest_distance > cur_city.costToGetTo(unchecked_cities[k]))
                        {
                            shortest_distance = cur_city.costToGetTo(unchecked_cities[k]);
                            index = k;
                        }

                    }

                    //This means no paths were found
                    if (shortest_distance == -1)
                    {
                        //I will do nothing and just go on to the next city in the list
                        break;
                    }
                    /*
                     * Then, once all the cities have been checked, I remove the shortest and check again until there's only one city left
                     */
                    //Set the current city I'm checking distances from and add it to the route
                    cur_city = unchecked_cities[index];
                    Route.Add(cur_city);

                    //Remove the current city from the list of cities to check
                    temp = new List<City>(unchecked_cities);
                    temp.RemoveAt(index);
                    unchecked_cities = temp.ToArray();

                    //reset the shortest distance and index variables
                    shortest_distance = -1;
                    index = 0;


                    //Add the last city at the end
                    if (unchecked_cities.Length == 1)
                    {
                        //If I can't get to the last city
                        if (!double.IsPositiveInfinity(cur_city.costToGetTo(unchecked_cities[0])))
                        {
                            //add the last city to the route
                            Route.Add(unchecked_cities[0]);
                            count++;
                            solution = new TSPSolution(Route);
                        }
                    }

                } while (unchecked_cities.Length > 1);

                //Console.Out.WriteLine("Cost of found route: " + solution.costOfRoute());
                //Console.Out.WriteLine("Cost of BSSF: " + costOfBssf());
                //Here I check the found route against the current best route
                if (solution != null)
                {
                    if (costOfBssf() == -1)
                    {
                        bssf = solution;
                    }
                    else if (costOfBssf() > solution.costOfRoute())
                    {
                        bssf = solution;
                    }
                }


            }
            timer.Stop();

            results[COST] = costOfBssf().ToString();
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = count.ToString();

            return results;
        }

        public string[] fancySolveProblem()
        {
            string[] results = new string[3];

            // TODO: Add your implementation for your advanced solver here.

            results[COST] = "not implemented";    // load results into array here, replacing these dummy values
            results[TIME] = "-1";
            results[COUNT] = "-1";

            return results;
        }
        #endregion

        public class matrix_and_bound
        {
            double[,] matrix;
            double lower_bound;

            public double[,] Matrix
            {
                get
                {
                    return matrix;
                }

                set
                {
                    matrix = value;
                }
            }

            public double Lower_bound
            {
                get
                {
                    return lower_bound;
                }

                set
                {
                    lower_bound = value;
                }
            }

            public matrix_and_bound(double[,] matrix, double lower_bound)
            {
                this.Matrix = matrix;
                this.Lower_bound = lower_bound;
            }
        }

        public class state_data
        {
            matrix_and_bound mb;
            ArrayList path = new ArrayList();
            ArrayList city_list = new ArrayList();
            int priority;
            int depth;
            int id;

            //These are constants I define to give weight to depth vs breadth
            int K = 3;
            int C = 5;

            public void add_city_index(int i)
            {
                city_list.Add(i);
            }
            public void add_city(City i)
            {
                path.Add(i);
            }

            public void set_priority()
            {
                priority = K * (int)mb.Lower_bound + C * depth;
            }

            public matrix_and_bound Mb
            {
                get
                {
                    return mb;
                }

                set
                {
                    mb = value;
                }
            }

            public ArrayList Path
            {
                get
                {
                    return path;
                }

                set
                {
                    path = value;
                }
            }

            public int Priority
            {
                get
                {
                    return priority;
                }

                set
                {
                    priority = value;
                }
            }

            public int Depth
            {
                get
                {
                    return depth;
                }

                set
                {
                    depth = value;
                }
            }

            public int Id
            {
                get
                {
                    return id;
                }

                set
                {
                    id = value;
                }
            }

            public ArrayList City_list
            {
                get
                {
                    return city_list;
                }

                set
                {
                    city_list = value;
                }
            }
        }

        public class PriorityQueue
        {
            state_data[] nodes; //the binary heap
            int[] pointers; //maps which node is at which index of nodes
            int next_open = 0;
            int length = 0;
            bool debug = false;

            public int Length
            {
                get
                {
                    return length;
                }

                set
                {
                    length = value;
                }
            }

            //Children are 2j + 1 and 2j + 2 (j is the index), the first element is the root and the last is the end
            public void make_queue(int num)
            {
                nodes = new state_data[num];
                pointers = new int[num];
                for (int i = 0; i < num; i++)
                {
                    nodes[i] = null;
                    pointers[i] = -1;
                }
            }
            public int get_num()
            {
                return next_open;
            }
            //just for testing purposes
            public state_data[] get_nodes()
            {
                return nodes;
            }
            public int[] get_pointers()
            {
                return pointers;
            }
            //Should run in O(log(|V|))
            public state_data delete_min()
            {
                if (debug) Console.WriteLine("\n DELETE MIN");
                if (debug) Console.WriteLine("\n Pointers[0] " + pointers[0]);
                //set the return index value
                state_data min = nodes[0];
                //change the first value to the last value and clear the last value
                pointers[0] = pointers[next_open - 1];
                nodes[0] = nodes[next_open - 1];
                if (debug) Console.WriteLine("nodes[next_open-1]: " + nodes[next_open - 1]);
                if (debug) Console.WriteLine("nodes[0] after switch: " + nodes[0]);

                nodes[next_open - 1] = null;
                pointers[next_open - 1] = -1;
                next_open--;
                if (debug) Console.WriteLine("next open: " + next_open);
                //sift_down the new root node
                sift_down(0);

                Length--;

                return min;
            }
            //Should run in O(log(|V|))
            public void insert(state_data data, int id)
            {
                if (debug) Console.WriteLine("INSERTING");
                pointers[next_open] = id;
                nodes[next_open] = data;
                bubble_up(next_open);
                next_open++;
                if (debug) Console.WriteLine("next open: " + next_open);
                Length++;
            }
            //This is what makes insertion take log time because the node might bubble up to the top of the tree
            private void bubble_up(int node)
            {
                if (debug) Console.WriteLine("BUBBLING");
                if (debug) Console.WriteLine("weight: " + nodes[node]);
                int p = ceiling(node);
                if (debug) Console.WriteLine("node: " + node + " p: " + p);
                if (debug) Console.WriteLine("node: " + nodes[node].Priority);
                if (debug && node != 0) Console.WriteLine("p: " + nodes[p].Priority);
                //While it's not the root and its parent is bigger than it is
                while (node != 0 && nodes[p].Priority > nodes[node].Priority)
                {
                    //now need to switch the nodes
                    state_data temp = nodes[node];
                    nodes[node] = nodes[p];
                    nodes[p] = temp;
                    //and update the pointers
                    int tempi = pointers[node];
                    pointers[node] = pointers[p];
                    pointers[p] = tempi;
                    //The current node is now the parent
                    node = p;
                    //update the parent
                    p = ceiling(node);
                    if (debug) Console.WriteLine("IN LOOP:: node: " + node + " p: " + p);
                }
                if (debug) Console.WriteLine("AT END:: node: " + node + " p: " + p);
                if (debug) Console.WriteLine("weight at end: " + nodes[node]);
            }
            //Gets the index of the node we want to change and updates the weight
            //The weight will always be a smaller value so we just change it and bubble up thenode, bubble up is what makes
             //it O(log(|V|))
            public void decrease_key(int node, int weight)
            {
                int k = 0;
                bool found = false;
                for (int i = 0; i < pointers.Length; i++)
                {
                    if (pointers[k] == node)
                    {
                        //we found the index in the array where node is stored
                        found = true;
                        break;
                    }
                    else
                    {
                        k++;
                    }
                }
                if (found) nodes[k].Priority = weight;
                            //if it wasn't found then it was already popped off the queue and doesn't needto be updated
                if (found) bubble_up(k);
            }

            //sift down will always work from the route (0) node
            public void sift_down(int node)
            {
                int c = get_min(node);
                if (debug) Console.WriteLine("SIFT DOWN get_min: " + c);
                if (c == -1)
                {
                    //it has no children
                }
                else
                {
                    while (c > 0 && nodes[c].Priority < nodes[node].Priority)
                    {
                        //now need to switch the nodes
                        state_data temp = nodes[node];
                        nodes[node] = nodes[c];
                        nodes[c] = temp;
                        //and update the pointers
                        int tempi = pointers[node];
                        pointers[node] = pointers[c];
                        pointers[c] = tempi;
                        //The current node is now the child
                        node = c;
                        //update the child
                        c = get_min(node);
                    }
                }
            }
            private int ceiling(int node)
            {
                double test = node;
                double parent = Math.Ceiling((test / 2) - 1);
                int p = (int)parent;
                return p;
            }
            public int get_min(int node)
            {
                //first check against the total number

                double test = node;
                int child_1 = (node * 2) + 1;
                int child_2 = (node * 2) + 2;
                //Make sure they're not out of bounds
                if (child_1 >= next_open)
                {
                    child_1 = -1;
                }
                if (child_2 >= next_open)
                {
                    child_2 = -1;
                }
                int c = child_1;
                if (child_1 == -1 && child_2 == -1)
                {
                    return c;
                }
                else if (child_1 == -1)
                {
                    return c;
                }
                else if (child_2 == -1)
                {
                    return c;
                }
                else if (nodes[child_1].Priority > nodes[child_2].Priority && child_2 != -1)
                {
                    c = child_2;
                }
                return c;
            }
        }

    }
}

    
