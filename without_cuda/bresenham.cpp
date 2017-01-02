
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace std;

vector<int> arange(int a, int b, int step)
{
  vector<int> ans;
  if(step == 1)
  {
    for(int i = a; i <= b; ++i)
    {
      ans.push_back(i);
    }
  }
  else
  {
    for(int i = a; i >= b; --i)
    {
      ans.push_back(i);
    }
  }
  return ans;
}

vector<int> cum(vector<int> q,  int x, char op)
{
  vector<int> ans;
  if(op == '+')
  {
    for(int i = 0; i < q.size(); ++i)
      {
        ans.push_back(q[i] + x );
      }
  }
  else
  {
    for(int i = 0; i < q.size(); ++i)
      {
        ans.push_back( x - q[i] );
      }
  }
  return ans;
}

void bresenham(  double x1,  double y1, double x2, double y2 )
{
  x1 = round(x1);
  x2 = round(x2);
  y1 = round(y1);
  y2 = round(y2);
  int dx = abs(x2 - x1);
  int dy = abs(y2 - y1);
  bool steep = abs(dy) > abs(dx);
  vector<int> q;

  if(steep)
  {
    int tmp = dx;
    dx = dy;
    dy = tmp;
  }
  if(dy == 0)
  {
    for(int  i = 0; i <= dx; ++i)
    {
      q.push_back(0);
    }
  }
  else
  {
    q.push_back(0);
    vector<int> arr;
    int cur;
    int cumSum = 0;
    for(int i = floor(dx/2); i >= floor(dx/2) - dy*dx; i -= dy)
    {
      cur = (i + dx)%dx;
      while(cur < 0)
      {
        cur = (cur + dx)%dx;
      }
      arr.push_back(cur);
     // cout<<cur<<endl;
    }
    for(int i = 1; i < arr.size(); ++i)
      {
        cur = arr[i] - arr[i-1];
        cur = cur >=0?1:0;
        cumSum += cur;
        q.push_back(cumSum);
       // cout<<cumSum<<endl;
      }
  }
  vector<int> y;
  vector<int> x;
  if(steep)
  {
    if(y1 <= y2)
    {
      y = arange(y1, y2, 1);
    }
    else
    {
      y = arange(y1, y2, -1);
    }

    if(x1 <= x2)
    {
      x = cum(q, x1, '+');
    }
    else
    {
      x = cum(q, x1, '-');
    }
  }
  else
  {
    if(x1 <= x2)
    {
      x = arange(x1, x2, 1);
    }
    else
    {
      x = arange(x1, x2, -1);
    }

    if(y1 <= y2)
    {
      y = cum(q, y1, '+');
    }
    else
    {
      y = cum(q, y1, '-');
    }
  }

  for(int i = 0; i < x.size(); ++i)
  {
        cout<<x[i]<<" "<<y[i]<<endl;
  }
  cout<<dx<< " " <<dy<<endl;
}


int main()
{
  double x1, x2, y1, y2;
  x1 = 1;
  x2 = 11;
  y1 = 2;
  y2 = 4;
  bresenham(x1, y1, x2, y2);
  return 0;
}