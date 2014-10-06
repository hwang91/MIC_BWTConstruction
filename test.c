//code1
/*
char* toStr() 
{
    char *s = "abcdefghijkl";
    return s;
}
int main()
{
    printf("%s\n",toStr());
    return 0;
}
*/
//code2
char* toStr() 
{
    char s[] = "abcdefghijkl";
    return s;
}
int main()
{
    printf("%s\n",toStr());
    return 0;
}
