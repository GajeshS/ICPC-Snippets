#include<bits/stdc++.h>
using namespace std;
#define ll long long
#define ld long double
#define mp make_pair
#define pb push_back
#define mod 1000000007
#define f first
#define s second
#define fastread ios_base::sync_with_stdio(false);cin.tie(NULL);cout.tie(NULL);
const int N=2e5+5;
const int maxn=1002;

ll segtree[4*N];
ll lazy[4*N];
ll cost[N];
/*
    @params
    cur: current node(0 indexed) left child is 2*cur+1 & right child is 2*cur+2
    lrange:left range of current node
    rrange:right range of current node
*/
void lazybuildSegTree(ll cur,ll lrange,ll rrange)
{
    if(lrange==rrange) //leaves
    {
        segtree[cur]=cost[lrange];
        return;
    }
    ll lc=2*cur+1,rc=2*cur+2;
    ll mid=(lrange+rrange)/2;
    lazybuildSegTree(lc,lrange,mid);
    lazybuildSegTree(rc,mid+1,rrange);
    segtree[cur]=segtree[lc]+segtree[rc]; // change here for different function
}

/*
    @params
    cur: current node number
    lrange,rrange: current node's range
    lq,rq: querying range
*/
ll lazyquerySegTree(ll cur,ll lrange,ll rrange,ll lq,ll rq)
{
    // Make pending updates using value stored in lazy 
    // nodes
    if(lazy[cur]!=0)
    {
        // Make pending updates using value stored in lazy 
        // nodes
        segtree[cur]+=((rrange-lrange+1)*lazy[cur]);   // change here for different function
        // checking if it is not leaf node because if 
        // it is leaf node then we cannot go further 
        if(rrange!=lrange) 
        {
            // We can postpone updating children we don't 
            // need their new values now. 
            // Since we are not yet updating children of si, 
            // we need to set lazy flags for the children 
            lazy[2*cur+1]+=lazy[cur];
            lazy[2*cur+2]+=lazy[cur];
        }
        lazy[cur]=0;
    }
    if(lrange>rrange)
        return 0;//non affecting value
    if(rrange<lq || lrange>rq)
        return 0;//non affecting value
    if(lrange>=lq && rrange<=rq)
        return segtree[cur];
    ll lc=2*cur+1,rc=2*cur+2;
    ll mid=(lrange+rrange)/2;
    ll q1=lazyquerySegTree(lc,lrange,mid,lq,rq);
    ll q2=lazyquerySegTree(rc,mid+1,rrange,lq,rq);
    return q1+q2; // change here for different function
}

/*
    range update
    @params
    cur: current node number
    lrange,rrange: current node's range
    ulrange,urrange: update query's range
    val: increment value
*/
void lazyrUpdate(ll cur,ll lrange,ll rrange,ll ulrange,ll urrange,ll val)
{
    // If lazy value is non-zero for current node of segment 
    // tree, then there are some pending updates. So we need 
    // to make sure that the pending updates are done before 
    // making new updates. Because this value may be used by 
    // parent after recursive calls (See last line of this 
    // function) 
    if(lazy[cur]!=0)
    {
        // Make pending updates using value stored in lazy 
        // nodes
        segtree[cur]+=((rrange-lrange+1)*lazy[cur]);   // change here for different function
        // checking if it is not leaf node because if 
        // it is leaf node then we cannot go further 
        if(rrange!=lrange) 
        {
            // We can postpone updating children we don't 
            // need their new values now. 
            // Since we are not yet updating children of si, 
            // we need to set lazy flags for the children 
            lazy[2*cur+1]+=lazy[cur];
            lazy[2*cur+2]+=lazy[cur];
        }
        lazy[cur]=0;
    }
    if(lrange>urrange || rrange<ulrange || lrange>rrange)
        return;
    if(lrange==rrange)
    {
        segtree[cur]+=val; // change here for different function
        return;
    }
    else if(lrange>=ulrange && rrange<=urrange) //completely includes
    {
        //incude the change in the current node of segtree
        segtree[cur]+=((rrange-lrange+1)*val);  // change here for different function
        /* postpone the update to the children */
        lazy[2*cur+1]+=val;
        lazy[2*cur+2]+=val;
        return;
    }
    ll lc=2*cur+1,rc=2*cur+2;
    ll mid=(lrange+rrange)/2;
    lazyrUpdate(lc,lrange,mid,ulrange,urrange,val);
    lazyrUpdate(rc,mid+1,rrange,ulrange,urrange,val);
    segtree[cur]=segtree[lc]+segtree[rc]; // change here for different function
}

void buildSegTree(ll cur,ll lrange,ll rrange)
{
    if(lrange==rrange) //leaves
    {
        segtree[cur]=cost[lrange];
        return;
    }
    ll lc=2*cur+1,rc=2*cur+2;
    ll mid=(lrange+rrange)/2;
    buildSegTree(lc,lrange,mid);
    buildSegTree(rc,mid+1,rrange);
    segtree[cur]=max(segtree[lc],segtree[rc]); // change here for different function
}

/*
    @params
    cur: current node number
    lrange,rrange: current node's range
    lq,rq: querying range
*/
ll querySegTree(ll cur,ll lrange,ll rrange,ll lq,ll rq)
{
    if(lrange>rrange)
        return 0;//non affecting value
    if(rrange<lq || lrange>rq)
        return 0;//non affecting value
    if(lrange>=lq && rrange<=rq)
        return segtree[cur];
    ll lc=2*cur+1,rc=2*cur+2;
    ll mid=(lrange+rrange)/2;
    ll q1=querySegTree(lc,lrange,mid,lq,rq);
    ll q2=querySegTree(rc,mid+1,rrange,lq,rq);
    return max(q1,q2); // change here for different function
}

/*
    point update
    @params
    cur: current node number
    lrange,rrange: current node's range
    pos: point update position
    val: update value(new value)
*/
void pUpdate(ll cur,ll lrange,ll rrange,ll pos,ll val)
{
    if(lrange>pos || rrange<pos || lrange>rrange)
        return;
    if(lrange==rrange)
    {
        segtree[cur]=val;
        return;
    }
    ll lc=2*cur+1,rc=2*cur+2;
    ll mid=(lrange+rrange)/2;
    pUpdate(lc,lrange,mid,pos,val);
    pUpdate(rc,mid+1,rrange,pos,val);
    segtree[cur]=max(segtree[lc],segtree[rc]); // change here for different function
}

//To compute x^y
//Complexity: log(n)
long long power(long long x,long long y)
{
    if(y==0)
        return 1;
    long long v=power(x,y/2);
    if(y%2==0)
    	return v*v;
    else
        return v*v*x;
}

ll gcd (ll n1, ll n2)
{
    if (n2 != 0)
       return gcd(n2, n1%n2);
    else 
       return n1;
}

/* Extended gcd algorithm
   Taken from: https://www.geeksforgeeks.org/euclidean-algorithms-basic-and-extended/
   solves ax+by=gcd(a,b)
*/
int gcdExtended(int a, int b, int *x, int *y) 
{ 
    // Base Case 
    if (a == 0) 
    { 
        *x = 0; 
        *y = 1; 
        return b; 
    } 
  
    int x1, y1; // To store results of recursive call 
    int gcd = gcdExtended(b%a, a, &x1, &y1); 
  
    // Update x and y using results of recursive 
    // call 
    *x = y1 - (b/a) * x1; 
    *y = x1; 
  
    return gcd; 
} 

//#define matrix_mod mod
template<class mat_type>
struct matrix {
 
    int n_rows, n_cols;
    vector< vector<mat_type> > m;
 
    matrix(int r = 1, int c = 1, bool I = false) {
        m.resize(r);
        for (int i = 0; i < r; i++) m[i].resize(c, 0);
        n_rows = r;
        n_cols = c;
        if (I)
            make_identity();
    }
 
    void make_identity() {
        for (int i = 0; i < n_rows; i++) m[i][i] = 1;
    }
 
#ifndef matrix_mod
 
    mat_type add(mat_type a, mat_type b) {
        return a + b;
    }
 
    mat_type sub(mat_type a, mat_type b) {
        return a - b;
    }
 
    mat_type mul(mat_type a, mat_type b) {
        return a * b;
    }
 
#else
 
    mat_type add(mat_type a, mat_type b) {
        return (a + b) % matrix_mod;
    }
 
    mat_type sub(mat_type a, mat_type b) {
        return ((a - b) % matrix_mod + matrix_mod) % matrix_mod;
    }
 
    mat_type mul(mat_type a, mat_type b) {
        return (a * b) % matrix_mod;
    }
 
#endif
 
    matrix operator +(const matrix& other) {
        matrix ans(n_rows, n_cols);
        for (int i = 0; i < n_rows; i++)
            for (int j = 0; j < n_cols; j++)
                ans.m[i][j] = add(m[i][j], other.m[i][j]);
        return ans;
    }
 
    matrix operator -(const matrix& other) {
        matrix ans(n_rows, n_cols);
        for (int i = 0; i < n_rows; i++)
            for (int j = 0; j < n_cols; j++)
                ans.m[i][j] = sub(m[i][j], other.m[i][j]);
        return ans;
    }
 
    matrix operator *(const matrix& other) {
        matrix ans(n_rows, other.n_cols);
        for (int i = 0; i < n_rows; i++)
            for (int j = 0; j < other.n_cols; j++)
                for (int k = 0; k < n_cols; k++)
                    ans.m[i][j] = add(ans.m[i][j], mul(m[i][k], other.m[k][j]));
        return ans;
    }
 
    matrix power(ll exp) {
        matrix ans(n_rows, n_cols, true);
        matrix multiplier = *this;
        while (exp > 0) {
            if (exp & 1)
                ans = ans * multiplier;
            exp >>= 1;
            multiplier = multiplier * multiplier;
        }
        return ans;
    }
	
};


ll add(ll a,ll b){
    return (a+b)%mod;
}
ll sub(ll a,ll b){
    return ((a-b)%mod+mod)%mod;
}
ll mul(ll a,ll b){
    return (a*b)%mod;
}

/* sieve of eratosthenes
   Taken from https://cp-algorithms.com/algebra/prime-sieve-linear.html
   Complexity: O(N log(logn))
   lp[i] stores the lowest prime number that is a divisor of i
   pr is the vector that stores the list of primes
*/
vector<ll>pr;
ll lp[N];
void sieve()
{
    for (int i=2; i<=N; ++i) {
        if (lp[i] == 0) {
            lp[i] = i;
            pr.push_back (i);
        }
        for (int j=0; j<(int)pr.size() && pr[j]<=lp[i] && i*pr[j]<=N; ++j)
            lp[i * pr[j]] = pr[j];
    }
}

bool compare(const pair<int,int>& lhs, const pair<int,int>& rhs)
{
  return lhs.first < rhs.first;
}

bool compare(const int* l,const int* r)
{
    return l<r;
}

ll fact[maxn],inv[maxn],ifact[maxn];
void calcInv(){
    for(int i=0;i<2;i++){
        fact[i]=1;
        ifact[i]=1;
        inv[i]=1;
    }
    int i;
    for(i=2;i<maxn;i++){
        fact[i]=(i*fact[i-1])%mod;
        ll q=mod/i;
        ll r=mod%i;
        inv[i]=((-q)*inv[r])%mod;
        inv[i]=(inv[i]+mod)%mod;
        ifact[i]=(inv[i]*ifact[i-1])%mod;
    }
}
		
ll ncr(ll n,ll r){
    if(n<r)
        return 0LL;
    else
    	return mul(mul(fact[n],ifact[r]),ifact[n-r]);
}

long long fermat_little(long long a,long long M)
{
    return power(a,M-2)%M;
}



void solve()
{
    return;
}

int main()
{ 
    fastread;
    int t;
    cin>>t;
    while(t--)
    {
        solve();
    }
    return 0;
}