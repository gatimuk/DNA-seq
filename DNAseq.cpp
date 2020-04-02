#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <string>
#include <cstdio>

using namespace std;

char codonToProtein(string codon)
{
  if (codon == "ATT" || codon == "ATC" || codon == "ATA") return 'I';
  if (codon == "CTT" || codon == "CTC" || codon == "CTA" || codon == "CTG" || codon == "TTA"|| codon == "TTG") return 'L';
  if (codon == "GTT" || codon == "GTC" || codon == "GTA" || codon == "GTG") return 'V';
  if (codon == "TTT" || codon == "TTC") return 'F';
  if (codon == "ATG") return 'M';
  if (codon == "TGT" || codon == "TGC") return 'C';
  if (codon == "GCT" || codon == "GCC" || codon == "GCA"|| codon == "GCG") return 'A';
  if (codon == "GGT" || codon == "GGC" || codon == "GGA"|| codon == "GGG") return 'G';
  if (codon == "CCT" || codon == "CCC" || codon == "CCA"|| codon == "CCG") return 'P';
  if (codon == "ACT" || codon == "ACC" || codon == "ACA"|| codon == "ACG") return 'T';
  if (codon == "TCT" || codon == "TCC" || codon == "TCA"|| codon == "TCG" || codon == "AGT"|| codon == "AGC") return 'S';
  if (codon == "TAT" || codon == "TAC") return 'Y';
  if (codon == "TGG") return 'W';
  if (codon == "CAA" || codon == "CAG") return 'Q';
  if (codon == "AAT" || codon == "AAC") return 'N';
  if (codon == "CAT" || codon == "CAC") return 'H';
  if (codon == "GAA" || codon == "GAG") return 'E';
  if (codon == "GAT" || codon == "GAC") return 'D';
  if (codon == "AAA" || codon == "AAG") return 'K';
  if (codon == "CGT" || codon == "CGC" || codon == "CGA" || codon == "CGG" || codon == "AGA" || codon == "AGG") return 'R';
  
  return 0;
}

bool isStartCodon(string codon)
{
  if (codon == "ATG") return true;
  return false;
}

bool isEndCodon(string codon)
{
  if (codon == "TAA" || codon == "TAG" || codon == "TGA") return true;
  return false;
}

////////////////////////

struct hashNode
{
  string key;
  int value;
  hashNode *next;
};

class hashTable
{
private:
  hashNode **hashes;
  unsigned int size;
public:
  hashTable(unsigned int _size) : size(_size) {
    hashes = new hashNode * [size];
  }

  hashTable() {
    hashes = new hashNode * [size];
  }

  ~hashTable() {
    delete[] hashes;
  }

  unsigned int hashFunc(string key) {
    unsigned long hash = 5381;
    for (size_t i = 0; i < key.size(); i++) {
      hash = ((hash << 5) + hash) + key[i];
    }
    return hash % size;
  } 

  void hashInsert(string key, int value){
    int index = hashFunc(key);
    hashNode* node = new hashNode;
    node->key = key;
    node->value = value;

    if (hashes[index] == NULL) {
      node->next = NULL;
      hashes[index] = node;
    } else {
      hashNode* temp = hashes[index];
      while (temp != NULL) {
        if (hashes[index]->key == temp->key) break;
        temp = temp->next;
      }
      if (temp == NULL) {
        node->next = hashes[index];
        hashes[index] = node;
      }
    }
  }

  bool isHashed(string key) {
    int index = hashFunc(key);
    if (hashes[index] == NULL) return false;
    else return true;
  }

  int getHashValue(string key) {
    int hashValue;
    int index = hashFunc(key);
    hashNode* node = hashes[index];
    while(node->next != NULL) node = node->next;
    hashValue = node->value;
    return hashValue;
  }
};

////////////////////////

struct seqObj
{
  string m_name;
  int m_freq;
};

struct geneObj : public hashTable
{
  string m_name;
  seqObj m_sequences[100];
  int numCodons;

  geneObj() : hashTable(1003) {
    numCodons = 0;
  }
};

////////////////////////


int main(int argc, char* argv[]) {
  //open file for input
  //ifstream inFile("input.txt");
  ifstream inFile(argv[1]);
  ofstream outFile("output.txt");
  char c;
  string totalSeq = "";
  string codonSeq = "";
  string geneSeq = "";
  bool isGene = false;

  int geneCounter = 0;
  int tableSize = 1003;
  geneObj geneRecords[100];
  hashTable ht(tableSize);

  //array of just genes
  string geneArray[101];

  unsigned int i;

  if (!inFile)
  {
    cerr << "unable to open input file\n";
    return -1;
  }

  while (inFile.get(c))
  {
    if (c == 'A' || c == 'C'|| c == 'G' || c == 'T' || c == 'a' || c == 'c'|| c == 'g' || c == 't')
    totalSeq += toupper(c);
  }

  //cout << totalSeq << endl;

  for(size_t i=0; i<totalSeq.size()-2; i++)
  {
    string codon = totalSeq.substr(i,3);
    char protein = codonToProtein(codon);

    if (!isGene && isStartCodon(codon))
    {
      isGene = true;
      codonSeq += codon;
      i += 2;
      continue;
    }

    //Single gene sequence constructed
    if (isGene && isEndCodon(codon))
    {
      isGene = false;
      geneSeq += protein;
      codonSeq += codon;

      if (ht.isHashed(geneSeq)) {
        //update pre-existing gene with codonSeq
        int geneIndex = ht.getHashValue(geneSeq);
        geneObj& gr = geneRecords[geneIndex];
        //update pre-exising codonSeq
        if (gr.isHashed(codonSeq))
        {
          int codonIndex = gr.getHashValue(codonSeq);
          gr.m_sequences[codonIndex].m_freq += 1;
        } else { //add new codon to pre-existing gene
          gr.hashInsert(codonSeq, gr.numCodons);
          gr.m_sequences[gr.numCodons].m_name = codonSeq;
          gr.m_sequences[gr.numCodons].m_freq += 1;
          gr.numCodons += 1;
        }
      } else {
        //add new gene
        ht.hashInsert(geneSeq, geneCounter);
        geneObj& gr = geneRecords[geneCounter];
        geneRecords[geneCounter].m_name = geneSeq;
        //add new codonSeq
        gr.hashInsert(codonSeq, gr.numCodons);
        gr.m_sequences[gr.numCodons].m_name = codonSeq;
        gr.m_sequences[gr.numCodons].m_freq += 1;
        gr.numCodons += 1;
        geneCounter++;  
      }
   
      geneSeq.erase();
      codonSeq.erase();
      i += 2;
      continue;
    }
    
    if (isGene)
    {
      geneSeq += protein;
      codonSeq += codon;
      i += 2;
    }
  }

  for(size_t i=0; i<geneCounter; i++) {
    //genObj& gr = geneRecords[i];
    outFile << geneRecords[i].m_name << endl ;
    for(size_t j=0; j<geneRecords[i].numCodons; j++) {
      outFile << "   " << geneRecords[i].m_sequences[j].m_name 
      << " " << geneRecords[i].m_sequences[j].m_freq << endl;
    }
  }

  inFile.close();
  outFile.close();
}