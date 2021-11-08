class Cell {
    private:
        // Methods
        void seed(long seed);
        long randl(long num);
        double randd(void);
        double gasdev(void);
        double variability(double mean);

        // Member variables
        double m_glucose;
        double m_dt;
        vector<double> m_results;

    public:
        // Constructor
        Cell(double glucose, double dt);
        
        // Destructor
        ~Cell();
        
        // Methods
        void set_initial(vector<double> initial);
        void update(void); // Dodati sklopitveni tok kot parameter
        void print_results(void);
        vector<double> get_result(void);
};