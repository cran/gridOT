#ifndef GRIDOT_TRANSPORT_STRUCT_H_INCLUDED
#define GRIDOT_TRANSPORT_STRUCT_H_INCLUDED

// encapsulates a transport from one index to another
struct Transport
{
	const int from;
	const int to;
	const double mass;

	Transport(int from, int to, double mass) : from(from), to(to), mass(mass)
	{
	}
};

#endif
