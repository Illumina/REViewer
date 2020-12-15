//
// REViewer
// Copyright 2020 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include "app/GenerateSvg.hh"

#include "graphcore/Path.hh"

#include <boost/optional.hpp>

#include <fstream>
#include <list>
#include <stdexcept>

using boost::optional;
using graphtools::NodeId;
using graphtools::Path;
using std::list;
using std::ofstream;
using std::string;
using std::to_string;
using std::vector;

static void
drawRect(ofstream& out, int x, int y, int width, int height, const string& fill, const string& stroke, double opacity)
{
    out << "<rect x=\"" << x << "\" y=\"" << y << "\"";
    out << " width=\"" << width << "\" height=\"" << height << "\"";
    out << " stroke=\"" << stroke << "\"";
    out << " fill=\"" << fill << "\"";
    out << " opacity=\"" << opacity << "\"";
    out << " />\n";
}

static void
drawRectWithLeftBreak(ofstream& out, int x, int y, int width, int height, const string& fill, const string& stroke)
{
    out << "<path d=\"";
    out << "M " << x << " " << y << " ";
    out << "h " << width << " ";
    out << "v " << height << " ";
    out << "h -" << width << " ";
    out << "l " << 5 << " " << -height / 4 << " ";
    out << "l " << -5 << " " << -height / 4 << " ";
    out << "l " << 5 << " " << -height / 4 << " ";
    out << "Z\"";

    out << " fill=\"" << fill << "\"";
    out << " stroke=\"" << stroke << "\"";

    out << "/>\n";
}

static void
drawRectWithRightBreak(ofstream& out, int x, int y, int width, int height, const string& fill, const string& stroke)
{
    out << "<path d=\"";
    out << "M " << x + width << " " << y << " ";
    out << "h " << -width << " ";
    out << "v " << height << " ";
    out << "h " << width << " ";
    out << "l " << -5 << " " << -height / 4 << " ";
    out << "l " << 5 << " " << -height / 4 << " ";
    out << "l " << -5 << " " << -height / 4 << " ";
    out << "Z\"";

    out << " fill=\"" << fill << "\"";
    out << " stroke=\"" << stroke << "\"";

    out << "/>\n";
}

static void drawLine(ofstream& out, int x, int y, int width, int height, const string& stroke)
{
    out << "<line ";
    out << "x1=\"" << x << "\" y1=\"" << y + height / 2 << "\" ";
    out << "x2=\"" << x + width << "\" y2=\"" << y + height / 2 << "\" ";
    out << "stroke=\"" << stroke << "\" ";
    out << "/>\n";
}

static void drawLetter(ofstream& out, int x, int y, int width, int height, char letter)
{
    string color = "black";
    if (letter == 'A')
    {
        color = "#FF6347";
    }
    else if (letter == 'T')
    {
        color = "#FCA100";
    }
    else if (letter == 'C')
    {
        color = "#393939";
    }
    else if (letter == 'G')
    {
        color = "#2F8734";
    }
    out << "<text x=\"" << x + width / 2 << "\" y=\"" << y + height / 2 << "\"";
    out << " dy=\"0.25em\""; // dominant-baseline="middle"
    out << " text-anchor=\"middle\"";
    out << " font-family=\"monospace\"";
    out << " font-size=\"11px\"";
    out << " stroke=\"" << color << "\"";
    out << ">";
    out << letter << "</text>";
}

static void drawText(ofstream& out, int x, int y, int width, int height, const string& text)
{
    const int letterWidth = width / text.length();
    for (int letterIndex = 0; letterIndex != text.length(); ++letterIndex)
    {
        drawLetter(out, x + letterWidth * letterIndex, y, letterWidth, height, text[letterIndex]);
    }
}

static void
drawArrows(ofstream& out, int x, int y, int width, int height, const string& stroke, const optional<string>& text)
{
    out << "<line x1=\"" << x << "\" y1=\"" << y + height / 2 << "\"";
    out << " x2=\"" << x + width << "\" y2=\"" << y + height / 2 << "\"";
    out << " stroke=\"" << stroke << "\"";
    out << " marker-start=\"url(#arrow)\" marker-end=\"url(#arrow)\" />\n";

    if (text)
    {
        out << "<text x=\"" << x + width / 2 << "\" y=\"" << y + height / 2 << "\"";
        out << " dy=\"0.25em\"";
        out << " text-anchor=\"middle\" font-family=\"monospace\" font-size=\"13px\"";
        out << " style=\"stroke:white; stroke-width:1.0em\" >";
        out << *text << "</text>\n";

        out << "<text x=\"" << x + width / 2 << "\" y=\"" << y + height / 2 << "\"";
        out << " dy=\"0.25em\"";
        out << " text-anchor=\"middle\" font-family=\"monospace\" font-size=\"13px\"";
        out << " font-weight=\"lighter\" stroke=\"" << stroke << "\" >";
        out << *text << "</text>\n";
    }
}

static void drawVerticalLine(ofstream& out, int x, int y, int width, int height, const string& stroke)
{
    out << "<line ";
    out << "x1=\"" << x << "\" y1=\"" << y << "\" ";
    out << "x2=\"" << x << "\" y2=\"" << y + height << "\" ";
    out << " style=\"stroke:" << stroke << "; stroke-width:3px\"";
    out << "/>\n";
}

void drawLane(ofstream& out, int baseWidth, int xPosStart, int yPos, const Lane& lane)
{
    for (const auto& segment : lane.segments)
    {
        int xPos = xPosStart + segment.start * baseWidth;
        for (const auto& feature : segment.features)
        {
            const int featureWidth = feature.length * baseWidth;
            if (feature.type == FeatureType::kRect)
            {
                drawRect(out, xPos, yPos, featureWidth, lane.height, feature.fill, feature.stroke, segment.opacity);
                if (feature.label)
                {
                    drawText(out, xPos, yPos, featureWidth, lane.height, *feature.label);
                }
            }
            else if (feature.type == FeatureType::kRectWithLeftBreak)
            {
                drawRectWithLeftBreak(out, xPos, yPos, featureWidth, lane.height, feature.fill, feature.stroke);
                if (feature.label)
                {
                    drawText(out, xPos, yPos, featureWidth, lane.height, *feature.label);
                }
            }
            else if (feature.type == FeatureType::kRectWithRightBreak)
            {
                drawRectWithRightBreak(out, xPos, yPos, featureWidth, lane.height, feature.fill, feature.stroke);
                if (feature.label)
                {
                    drawText(out, xPos, yPos, featureWidth, lane.height, *feature.label);
                }
            }
            else if (feature.type == FeatureType::kLine)
            {
                drawLine(out, xPos, yPos, featureWidth, lane.height, feature.stroke);
            }
            else if (feature.type == FeatureType::kArrows)
            {
                drawArrows(out, xPos, yPos, featureWidth, lane.height, feature.stroke, feature.label);
            }
            else if (feature.type == FeatureType::kVerticalLine)
            {
                drawVerticalLine(out, xPos, yPos, featureWidth, lane.height, feature.stroke);
            }
            else
            {
                throw std::runtime_error("Encountered feature of unknown type");
            }

            xPos += featureWidth;
        }
    }
}

int getPlotWidth(int baseWidth, const vector<LanePlot>& lanePlots)
{
    int maxLaneWidth = 0;
    for (const auto& lanePlot : lanePlots)
    {
        for (const auto& lane : lanePlot)
        {
            for (const auto& segment : lane.segments)
            {
                maxLaneWidth = std::max(maxLaneWidth, segment.end);
            }
        }
    }

    return maxLaneWidth * baseWidth;
}

int getPlotHeight(int spacingBetweenLanes, int spacingBetweenLanePlots, const vector<LanePlot>& lanePlots)
{
    int plotHeight = 0;
    for (const auto& lanePlot : lanePlots)
    {
        if (plotHeight != 0)
        {
            plotHeight += spacingBetweenLanePlots;
        }
        for (const auto& lane : lanePlot)
        {
            plotHeight += lane.height;
        }
        plotHeight += static_cast<int>(lanePlot.size()) * spacingBetweenLanes;
    }

    return plotHeight;
}

void generateSvg(const vector<LanePlot>& lanePlots, const string& outputPath)
{
    const int kSpacingBetweenLanes = 5;
    const int kSpacingBetweenLanePlots = 50;
    const int kBaseWidth = 10;
    const int kPlotPadX = 10;
    const int kPlotPadY = 5;

    const int plotWidth = getPlotWidth(kBaseWidth, lanePlots) + 2 * kPlotPadX;
    const int plotHeight = getPlotHeight(kSpacingBetweenLanes, kSpacingBetweenLanePlots, lanePlots) + 2 * kPlotPadY;

    int yPos = kPlotPadY;

    ofstream svgFile(outputPath);
    if (svgFile.is_open())
    {
        svgFile << "<svg width=\"" << plotWidth << "\" height=\"" << plotHeight << "\""
                << " xmlns=\"http://www.w3.org/2000/svg\">\n";
        svgFile << "<defs>\n"
                   "    <linearGradient id=\"BlueWhiteBlue\" x1=\"0%\" y1=\"0%\" x2=\"0%\" y2=\"100%\">\n"
                   "      <stop offset=\"0%\" style=\"stop-color:#8da0cb;stop-opacity:0.8\" />\n"
                   "      <stop offset=\"50%\" style=\"stop-color:#8da0cb;stop-opacity:0.1\" />\n"
                   "      <stop offset=\"100%\" style=\"stop-color:#8da0cb;stop-opacity:0.8\" />\n"
                   "    </linearGradient>\n"
                   "    <linearGradient id=\"OrangeWhiteOrange\" x1=\"0%\" y1=\"0%\" x2=\"0%\" y2=\"100%\">\n"
                   "      <stop offset=\"0%\" style=\"stop-color:#fc8d62;stop-opacity:0.8\" />\n"
                   "      <stop offset=\"50%\" style=\"stop-color:#fc8d62;stop-opacity:0.1\" />\n"
                   "      <stop offset=\"100%\" style=\"stop-color:#fc8d62;stop-opacity:0.8\" />\n"
                   "    </linearGradient>\n"
                   "    <linearGradient id=\"GreenWhiteGreen\" x1=\"0%\" y1=\"0%\" x2=\"0%\" y2=\"100%\">\n"
                   "      <stop offset=\"0%\" style=\"stop-color:#66c2a5;stop-opacity:0.8\" />\n"
                   "      <stop offset=\"50%\" style=\"stop-color:#66c2a5;stop-opacity:0.1\" />\n"
                   "      <stop offset=\"100%\" style=\"stop-color:#66c2a5;stop-opacity:0.8\" />\n"
                   "    </linearGradient>\n"
                   "    <marker id=\"arrow\" viewBox=\"0 0 10 10\" refX=\"9\" refY=\"5\"\n"
                   "        markerWidth=\"6\" markerHeight=\"6\"\n"
                   "        orient=\"auto-start-reverse\">\n"
                   "      <path d=\"M 0 0 L 10 5 L 0 10 z\" />\n"
                   "    </marker>"
                   "</defs>";

        for (const auto& lanePlot : lanePlots)
        {
            for (const auto& lane : lanePlot)
            {
                drawLane(svgFile, kBaseWidth, kPlotPadX, yPos, lane);
                yPos += lane.height + kSpacingBetweenLanes;
            }

            yPos += kSpacingBetweenLanePlots;
        }

        svgFile << "</svg>" << std::endl;
    }
    else
    {
        throw std::runtime_error("Unable to open " + outputPath);
    }
}
