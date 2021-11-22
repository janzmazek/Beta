class Cell {
    private:
        // Methods
        void seed(long seed);
        long randl(long num);
        double randd(void);
        double gasdev(void);
        double variability(double mean, double std, mt19937 generator);

        // Static member variables
        double m_glucose, m_dt;
        vector<double> m_results;

        // Variable member variables
        double m_cm, m_gca, m_gkatp, m_gk, m_gs, m_gkatp50, m_alfa, m_kca;

    public:
        // Constructor
        Cell(double glucose, double dt, double initial[6], mt19937 generator);
        
        // Destructor
        ~Cell();
        
        // Methods
        void update(double coupling_parameter); // Dodati sklopitveni tok kot parameter
        void print_results(void);
        vector<double> get_result(void);
};