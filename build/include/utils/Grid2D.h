#ifndef GRID2D_H
#define GRID2D_H

#include <map>
#include "TMath.h"

namespace Utils {
    template <typename T>
    class Grid2D {
    private:
        std::map<std::pair<UInt_t, UInt_t>, T> grid;
        UInt_t nrows;
        UInt_t ncols;
    public:
        T &at(UInt_t i, UInt_t j) {
            return grid[std::make_pair(i, j)];
        }

        std::map<UInt_t, T> row(UInt_t i) {
            std::map<UInt_t, T> result;
            for (const auto& kv : grid) {
                if (kv.first.first == i) {
                    result[kv.first.second] = kv.second;
                }
            }
            return result;
        }

        std::map<UInt_t, T> col(UInt_t j) {
            std::map<UInt_t, T> result;
            for (const auto& kv : grid) {
                if (kv.first.second == j) {
                    result[kv.first.first] = kv.second;
                }
            }
            return result;
        }

        std::vector<UInt_t> jsInRow(UInt_t i) {
            std::vector<UInt_t> result;
            for (const auto& kv : grid) {
                if (kv.first.first == i) {
                    result.push_back(kv.first.second);
                }
            }
            return result;
        }

        std::vector<UInt_t> isInCol(UInt_t j) {
            std::vector<UInt_t> result;
            for (const auto& kv : grid) {
                if (kv.first.second == j) {
                    result.push_back(kv.first.first);
                }
            }
            return result;
        }

        std::vector<UInt_t> rowIs() {
            std::vector<UInt_t> result;
            for (const auto& kv : grid) {
                if (std::find(result.begin(), result.end(), kv.first.first) == result.end()) {
                    result.push_back(kv.first.first);
                }
            }
            return result;
        }

        std::vector<UInt_t> colJs() {
            std::vector<UInt_t> result;
            for (const auto& kv : grid) {
                if (std::find(result.begin(), result.end(), kv.first.second) == result.end()) {
                    result.push_back(kv.first.second);
                }
            }
            return result;
        }

        void clearRowExceptOne(UInt_t rowI) {
            std::map<UInt_t, Float_t> row = this->row(rowI);
            std::map<UInt_t, Float_t>::iterator it = row.begin();
            ++it;
            for ( ; it != row.end(); ++it) {
                grid.erase(std::make_pair(rowI, (*it).first));
            }
        }

        void clearColExceptOne(UInt_t colJ) {
            std::map<UInt_t, Float_t> col = this->col(colJ);
            std::map<UInt_t, Float_t>::iterator it = col.begin();
            ++it;
            for ( ; it != col.end(); ++it) {
                grid.erase(std::make_pair((*it).first, colJ));
            }
        }

        bool hasI(UInt_t i) {
            for (const auto& kv : grid) {
                if (kv.first.first == i) {
                    return true;
                }
            }
            return false;
        }

        bool hasJ(UInt_t j) {
            for (const auto& kv : grid) {
                if (kv.first.second == j) {
                    return true;
                }
            }
            return false;
        }

        bool canReshuffleToDiagonal() {
            std::vector<UInt_t> rowIs = rowIs();
            for (UInt_t rowI : rowIs) {
                std::map<UInt_t, T> row = row(rowI);
                if (row.size() != 1) {
                    return false;
                }
            }
            return true;
        }
    };
}
#endif
