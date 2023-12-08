#include <cstdlib>
#include <iostream>
#include <mpi.h>
#include <time.h>

using namespace std;

class Position {
public:
  int line;
  int column;
  Position(int line, int column) : line(line), column(column) {}
};

class Cell {
public:
  int N;
  int line;
  int column;

  Cell(int N, int line, int column) : N(N), line(line), column(column) {}

  Position getLeft() {
    int _column = column == 0 ? (N - 1) : (column - 1);
    return Position(this->line, _column);
  }

  Position getRight() {
    int _column = column == (N - 1) ? 0 : (column + 1);
    return Position(this->line, _column);
  }

  Position getUpper() {
    int _line = line == 0 ? (N - 1) : (line - 1);
    return Position(_line, this->column);
  }

  Position getLower() {
    int _line = line == (N - 1) ? 0 : (line + 1);
    return Position(_line, this->column);
  }
};

class Grid {
public:
  float **grid;
  int N;

  Grid(int N) : N(N) {
    grid = new float *[N];
    for (int lineIndex = 0; lineIndex < N; lineIndex++) {
      grid[lineIndex] = new float[N];
      for (int columnIndex = 0; columnIndex < N; columnIndex++) {
        grid[lineIndex][columnIndex] = 0;
      }
    }
  }

  ~Grid() {
    for (int lineIndex = 0; lineIndex < N; lineIndex++) {
      delete[] grid[lineIndex];
    }
    delete[] grid;
  }

  float getCellValue(int line, int column) { return this->grid[line][column]; }

  Cell getCell(int line, int column) { return Cell(N, line, column); }

  float getCellNeighboursAverage(int line, int column) {
    Cell cell = this->getCell(line, column);
    Position upper = cell.getUpper();
    Position lower = cell.getLower();
    Position left = cell.getLeft();
    Position right = cell.getRight();
    Position positions[] = {upper,
                            lower,
                            left,
                            right,
                            Position(upper.line, left.column),
                            Position(upper.line, right.column),
                            Position(lower.line, left.column),
                            Position(lower.line, right.column)};
    int sum = 0;
    for (int index = 0; index < 8; index++) {
      Position position = positions[index];
      sum += this->getCellValue(position.line, position.column);
    }
    return sum / 8.0;
  }

  int countCellNeighbours(int line, int column) {
    Cell cell = this->getCell(line, column);
    Position upper = cell.getUpper();
    Position lower = cell.getLower();
    Position left = cell.getLeft();
    Position right = cell.getRight();
    Position positions[] = {upper,
                            lower,
                            left,
                            right,
                            Position(upper.line, left.column),
                            Position(upper.line, right.column),
                            Position(lower.line, left.column),
                            Position(lower.line, right.column)};
    int count = 0;
    for (int index = 0; index < 8; index++) {
      Position position = positions[index];
      if (this->getCellValue(position.line, position.column))
        count++;
    }
    return count;
  }

  float getNextCellState(int line, int column) {
    int neighbours = this->countCellNeighbours(line, column);
    float cellValue = this->getCellValue(line, column);
    if (cellValue > 0 && (neighbours == 2 || neighbours == 3))
      return cellValue;
    if (cellValue == 0 && neighbours == 3)
      return 1;
    return 0;
  }

  Grid *getNewGrid(int startLine, int endLine) {
    Grid *newGrid = new Grid(N);
    for (int lineIndex = startLine; lineIndex < endLine; lineIndex++) {
      for (int columnIndex = 0; columnIndex < N; columnIndex++) {
        newGrid->grid[lineIndex][columnIndex] =
            this->getNextCellState(lineIndex, columnIndex);
      }
    }
    return newGrid;
  }

  void updateBorders(int startLine, int endLine, int processRank, Grid *grid) {
    MPI_Status status;

    // Send the bottom row to the process below
    if (endLine < grid->N) {
      MPI_Sendrecv(grid->grid[endLine - 1], grid->N, MPI_FLOAT, processRank + 1,
                   0, grid->grid[endLine], grid->N, MPI_FLOAT, processRank + 1,
                   0, MPI_COMM_WORLD, &status);
    }

    // Send the top row to the process above
    if (startLine > 0) {
      MPI_Sendrecv(grid->grid[startLine], grid->N, MPI_FLOAT, processRank - 1,
                   0, grid->grid[startLine - 1], grid->N, MPI_FLOAT,
                   processRank - 1, 0, MPI_COMM_WORLD, &status);
    }
  }

  int countRemainingLivingCells() {
    int livingCells = 0;
    for (int l = 0; l < N; l++) {
      for (int c = 0; c < N; c++) {
        if (this->grid[l][c] > 0)
          livingCells++;
      }
    }
    return livingCells;
  }

  void printAll() {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        if (this->grid[i][j] == 0) {
          cout << "□";
        } else {
          cout << "■";
        }
      }
      cout << "\n";
    }
  }

  void print() {
    int max = N >= 50 ? 50 : N;
    cout << "  ";
    for (int i = 0; i < max; i++) {
      if (i % 2 == 0) {
        cout << "▬▬";
      } else {
        cout << "▭▭";
      }
    }
    cout << "\n";
    for (int i = 0; i < max; i++) {
      if (i % 2 == 0) {
        cout << " ▌";
      } else {
        cout << "  ";
      }
      for (int j = 0; j < max; j++) {
        if (this->grid[i][j] == 0) {
          cout << "  ";
        } else {
          cout << "■■";
        }
      }
      if (i % 2 == 0) {
        cout << "▌";
      } else {
        cout << "  ";
      }
      cout << "\n";
    }
    cout << "  ";
    for (int i = 0; i < max; i++) {
      if (i % 2 == 0) {
        cout << "▬▬";
      } else {
        cout << "▭▭";
      }
    }
    cout << "\n";
  }
};

Grid *play(int generations, int numProcesses, int processRank, Grid *grid) {
  int i;
  Grid *newGrid, *_grid = grid->getNewGrid(0, grid->N);
  cout << "Generation: 1 - " << _grid->countRemainingLivingCells() << endl;

  for (i = 0; i < generations - 1; i++) {
    // Update borders before calculating the new grid
    _grid->updateBorders(0, _grid->N, processRank, _grid);
    newGrid = _grid->getNewGrid(0, grid->N);
    delete _grid;
    _grid = newGrid;
    cout << "Generation: " << i + 2 << " - "
         << _grid->countRemainingLivingCells() << " living cells." << endl;
  }

  return _grid;
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  double start, end;
  int numProcesses, processRank;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
  MPI_Comm_rank(MPI_COMM_WORLD, &processRank);

  start = MPI_Wtime();

  int generations = 2000, N = 2048;
  int linesPerProcess = N / numProcesses;

  Grid *grid = new Grid(N);

  // Initialize the table in parallel
  srand(time(NULL) + processRank);
  for (int i = processRank * linesPerProcess;
       i < (processRank + 1) * linesPerProcess; i++) {
    for (int j = 0; j < N; j++) {
      grid->grid[i][j] = rand() % 2;
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  Grid *finalGrid = play(generations, numProcesses, processRank, grid);

  delete finalGrid;

  MPI_Barrier(MPI_COMM_WORLD);


  // print results
  if (processRank == 0) {

    cout << "\nThere are still " << grid->countRemainingLivingCells()
         << " living cells." << endl;
    end = MPI_Wtime();
    cout << "Took " << end - start << " seconds." << endl;
  }

  MPI_Finalize();
  return 0;
}
