#include "lbm.h"

void gaussian_density(int center, int std, arma::mat &rho)
{
    for (arma::uword x = 0; x < WIDTH; x++)
    {
        for (arma::uword y = 0; y < HEIGHT; y++)
        {

            double modifier = std::pow(2.7, -1.0f * std::pow(double(x) - center, 2.0f) / double(std)) * std::pow(2.7, -1.0f * std::pow(double(y) - center, 2.0f) / float(std));
            rho(x, y) += modifier;
        }
    }
}

int main()
{
    std::cout << omega << std::endl;
    std::cout << Re << std::endl;
    std::cout << C * dt / dx << std::endl;
    arma::cube fin = arma::cube(WIDTH, HEIGHT, Q);
    arma::cube fout = arma::cube(WIDTH, HEIGHT, Q);
    arma::cube feq = arma::cube(WIDTH, HEIGHT, Q);

    arma::mat rho = arma::mat(WIDTH, HEIGHT);
    arma::cube u_initial = arma::cube(WIDTH, HEIGHT, 2); // initial velocity field
    arma::cube u = arma::cube(WIDTH, HEIGHT, 2);

    // init velocity or density
    u_initial.zeros();
    //u_initial.randu();
    //u_initial -= .5f;
    //u_initial *= ulb;
    u = u_initial;

    rho.ones();
    gaussian_density(WIDTH / 2, WIDTH / 10, rho);
    rho /= rho.max();

    equilibrium(rho, u, feq);
    fout = feq;

    for (int i = 0; i < MAX_TIMES * steps; i++) // iterate over given number of characteristic times
    {
        if (!(i % TIME))
        {
            std::ostringstream name;
            name << "frame" << i / TIME << ".txt";
            u.save(name.str(), arma::arma_ascii);
        }

        stream(fout, fin); // naiively advances all states

        // find density and velocity
        sumpop(fin, rho);
        velocity(fin, rho, u);

        equilibrium(rho, u, feq);
        collide(fin, feq, fout);
    }

    return 0;
}