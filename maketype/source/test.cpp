#include <vector>
#include <iostream>

using namespace std;

void func(vector<int> a)
{
  cout << "good" << endl;
}

int main(int argc, char** argv)
{
  vector< vector<int> > a = {{1,2},{3,4}};
  func(a[0]);
}
