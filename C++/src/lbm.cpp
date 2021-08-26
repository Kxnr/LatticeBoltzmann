#include "lbm.h"

LatticeSimulator LatticeSimulator step() {
    return this.stream().bounceback().equilibrate()
}

D2Q9LatticeSimulator LatticeSimulator equilibrate() {
    // calculate equilibrium distribution

    Lattice density =

    // sum over (distribution * velocity) / density

    eq_lattice = // TODO

    // run collide step with relaxation towards equilibrium
    return this.copy().set_lattice(this.lattice - this.omega * (this.lattice - eq_lattice));
}

void equilibrium(arma::mat rho, arma::cube u, arma::cube &feq)
{
    arma::cube cv(WIDTH, HEIGHT, Q);
    for (arma::uword q = 0; q < Q; q++)
    {
        cv.slice(q) = 3.0f * (u.slice(0) * vel(0, q) + u.slice(1) * vel(1, q));
    }

    arma::mat usqr = arma::sum(arma::pow(u, 2), 2);

    for (arma::uword q = 0; q < Q; q++)
    {
        feq.slice(q) = rho % (w(q) * (1.0f + cv.slice(q) + .5f * arma::pow(cv.slice(q), 2) - (3.0f / 2.0f) * usqr));
    }
}

void sumpop(arma::cube pop, arma::mat &rho)
{
    rho.zeros();
    for (arma::uword q = 0; q < Q; q++)
    {
        rho += pop.slice(q);
    }
}

template <size_t S>
void bounceback(arma::cube fin, arma::uword (*indices)[S][2], arma::cube &fout)
{
    for (arma::uword q = 0; q < Q; q++)
    {
        for (int i = 0; i < S; i++)
        {
            fout(indices[0], indices[1], q) = fin(indices[0], indices[1], opp[q]);
        }
    }
}

void collide(arma::cube fin, arma::cube feq, arma::cube &fout)
{
    fout = fin - omega * (fin - feq);
}

void stream(arma::cube fin, arma::cube &fout)
{
    arma::cube temp(WIDTH, HEIGHT, Q);
    arma::sword xoff, yoff, idx, idy;
    for (arma::uword q = 0; q < Q; q++)
    {
        xoff = vel(0, q);
        for (arma::sword x = 0; x < WIDTH; x++)
        {
            idx = (x + xoff) % WIDTH;
            idx = (idx < 0) ? WIDTH - 1 : idx;

            temp.slice(q).row(idx) = fin.slice(q).row(x);
        }
    }
    for (arma::uword q = 0; q < Q; q++)
    {
        yoff = vel(1, q);
        for (arma::sword y = 0; y < HEIGHT; y++)
        {
            idy = (y + yoff) % HEIGHT;
            idy = (idy < 0) ? HEIGHT - 1 : idy;

            fout.slice(q).col(idy) = temp.slice(q).col(y);
        }
    }
}

void velocity(arma::cube fin, arma::mat rho, arma::cube &u)
{
    arma::cube temp(WIDTH, HEIGHT, Q);
    for (arma::uword c = 0; c < 2; c++)
    {
        for (arma::uword q = 0; q < Q; q++)
        {
            temp.slice(q) = fin.slice(q) * vel(c, q);
        }
        u.slice(c) = arma::sum(temp, 2);
        u.slice(c) /= rho;
    }
}

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
