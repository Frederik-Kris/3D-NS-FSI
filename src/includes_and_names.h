/*
 * includes_and_names.h
 *
 *  Created on: Mar 15, 2021
 *      Author: frederik
 *
 *  The purpose of this file is to have all the library inclusions and 'using' statements in one place,
 *  instead of writing 'using std::cout;', etc., in all the header files.
 */

#ifndef SRC_INCLUDES_AND_NAMES_H_
#define SRC_INCLUDES_AND_NAMES_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <set>
#include <exception>
#include <math.h>
#include <algorithm>
#include <SFML/Graphics.hpp>
#include <libconfig.h++>


using std::cout;
using std::ofstream;
using std::endl;
using std::setprecision;
using std::string;
using std::to_string;
using std::vector;
using std::unique_ptr;
using std::invalid_argument;
using std::pow;
using std::fabs;
using std::sqrt;
using std::max;
using std::min;
using std::fmod;

using sf::RenderWindow;
using sf::VideoMode;
using sf::Event;
using sf::Keyboard;
using sf::RectangleShape;
using sf::Vector2f;
using sf::Color;
using sf::Clock;
using sf::Time;

using libconfig::Config;
using libconfig::Setting;


typedef std::array<size_t, 8> Array8_u;
typedef std::array<bool,   8> Array8_b;


#endif /* SRC_INCLUDES_AND_NAMES_H_ */
