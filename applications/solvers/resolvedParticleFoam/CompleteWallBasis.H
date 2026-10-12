// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef LPBF_COMPLETE_WALL_BASIS_H
#define LPBF_COMPLETE_WALL_BASIS_H
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
namespace lpbfM2
{
struct CompleteWallBasis
{
    std::vector<std::vector<double>> rows,coordinates;
    double coverage=0,orthogonality=0;
    explicit CompleteWallBasis(const std::vector<std::vector<double>>& original)
    {
        if(original.empty()||original[0].empty()) throw std::runtime_error("Empty wall sampling space");
        const std::size_t columns=original[0].size();
        auto residual=original;
        auto dot=[](const std::vector<double>& a,const std::vector<double>& b)
        {double value=0;for(std::size_t i=0;i<a.size();++i) value+=a[i]*b[i];return value;};
        for(const auto& row:original)
        {if(row.size()!=columns) throw std::runtime_error("Wall sampling dimensions differ");
         for(double v:row) if(!std::isfinite(v)) throw std::runtime_error("Nonfinite wall sampling");}
        for(std::size_t step=0;step<original.size();++step)
        {
            std::size_t pivot=0;double maximum=-1;
            for(std::size_t i=0;i<residual.size();++i)
            {double norm=dot(residual[i],residual[i]);if(norm>maximum){maximum=norm;pivot=i;}}
            auto value=residual[pivot];
            for(int pass=0;pass<2;++pass) for(const auto& row:rows)
            {double projection=dot(value,row);for(std::size_t k=0;k<columns;++k) value[k]-=projection*row[k];}
            double norm=std::sqrt(dot(value,value));
            if(norm<=1e-11) break;
            for(double& v:value) v/=norm;
            rows.push_back(value);
            for(auto& row:residual) for(int pass=0;pass<2;++pass)
            {double projection=dot(row,value);for(std::size_t k=0;k<columns;++k) row[k]-=projection*value[k];}
        }
        if(rows.empty()) throw std::runtime_error("Zero wall sampling rank");
        for(const auto& row:original)
        {
            std::vector<double> coordinate,remainder=row;
            for(const auto& basis:rows)
            {double coefficient=dot(row,basis);coordinate.push_back(coefficient);
             for(std::size_t k=0;k<columns;++k) remainder[k]-=coefficient*basis[k];}
            coverage=std::max(coverage,std::sqrt(dot(remainder,remainder)/dot(row,row)));
            coordinates.push_back(coordinate);
        }
        for(std::size_t i=0;i<rows.size();++i) for(std::size_t j=0;j<=i;++j)
            orthogonality=std::max(orthogonality,std::abs(dot(rows[i],rows[j])-(i==j ? 1.:0.)));
        if(!std::isfinite(coverage)||coverage>1e-10||orthogonality>1e-10)
            throw std::runtime_error("Complete wall basis coverage/orthogonality failed");
    }
};
}
#endif
