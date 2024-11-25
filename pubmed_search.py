import streamlit as st
from Bio import Entrez
from typing import List, Dict
import pandas as pd
from datetime import datetime
import time

class PubMedSearch:
    def __init__(self, email: str, api_key: str):
        """
        Inicializa o buscador com credenciais
        """
        self.email = email
        self.api_key = api_key
        Entrez.email = email
        Entrez.api_key = api_key

    def search_articles(self,
                        query: str,
                        max_results: int = 10,
                        start_year: int = None,
                        end_year: int = None) -> List[Dict]:
        """
        Realiza busca no PubMed
        """
        # Adiciona filtro de data se especificado
        if start_year and end_year:
            query = f"{query} AND ({start_year}:{end_year}[pdat])"

        try:
            with st.spinner('Buscando artigos...'):
                # Busca IDs dos artigos
                search_handle = Entrez.esearch(
                    db="pubmed",
                    term=query,
                    retmax=max_results,
                    sort="relevance"
                )
                search_results = Entrez.read(search_handle)
                search_handle.close()

                if not search_results["IdList"]:
                    return []

                # Busca detalhes dos artigos
                fetch_handle = Entrez.efetch(
                    db="pubmed",
                    id=search_results["IdList"],
                    rettype="abstract",
                    retmode="xml"
                )

                articles = Entrez.read(fetch_handle)['PubmedArticle']
                fetch_handle.close()

                # Processa resultados
                processed_articles = []
                for article in articles:
                    try:
                        article_data = article['MedlineCitation']['Article']

                        # Extrai autores
                        authors = []
                        if 'AuthorList' in article_data:
                            for author in article_data['AuthorList']:
                                if 'LastName' in author and 'ForeName' in author:
                                    authors.append(f"{author['LastName']} {author['ForeName']}")

                        # Extrai abstract
                        abstract = ""
                        if 'Abstract' in article_data:
                            abstract = article_data['Abstract']['AbstractText'][0]

                        # Extrai DOI
                        doi = ""
                        if 'ELocationID' in article_data:
                            for id in article_data['ELocationID']:
                                if id.attributes['EIdType'] == 'doi':
                                    doi = str(id)

                        processed_articles.append({
                            'title': article_data['ArticleTitle'],
                            'authors': '; '.join(authors[:3]) + (' et al.' if len(authors) > 3 else ''),
                            'journal': article_data['Journal']['Title'],
                            'year': int(article_data['Journal']['JournalIssue']['PubDate'].get('Year', 0)),
                            'abstract': abstract[:300] + '...' if len(abstract) > 300 else abstract,
                            'pmid': article['MedlineCitation']['PMID'].contents[0],
                            'doi': doi,
                            'url': f"https://pubmed.ncbi.nlm.nih.gov/{article['MedlineCitation']['PMID'].contents[0]}"
                        })

                    except Exception as e:
                        st.error(f"Erro ao processar artigo: {str(e)}")
                        continue

                return processed_articles

        except Exception as e:
            st.error(f"Erro na busca: {str(e)}")
            return []

def main():
    st.title("🔍 Buscador PubMed")
    st.write("Pesquise artigos científicos na base do PubMed")

    # Configurações na sidebar
    with st.sidebar:
        st.header("Configurações")
        email = st.text_input("Email", key="email")
        api_key = st.text_input("API Token PubMed", type="password", key="api_key")

        st.markdown("---")
        st.markdown("""
        ### Dicas de Pesquisa
        - Use AND para combinar termos
        - Use OR para alternativas
        - Use aspas para frases exatas
        - Use [MeSH Terms] para termos MeSH
        """)

    # Formulário de busca
    with st.form("search_form"):
        col1, col2 = st.columns([3, 1])

        with col1:
            query = st.text_input("Termo de busca",
                                  placeholder="Ex: diabetes AND nutrition therapy[MeSH Terms]")

        with col2:
            max_results = st.number_input("Máx. resultados",
                                          min_value=1, max_value=100, value=10)

        col3, col4 = st.columns(2)
        with col3:
            start_year = st.number_input("Ano inicial",
                                         min_value=1900,
                                         max_value=datetime.now().year,
                                         value=2020)
        with col4:
            end_year = st.number_input("Ano final",
                                       min_value=1900,
                                       max_value=datetime.now().year,
                                       value=datetime.now().year)

        submitted = st.form_submit_button("Buscar")

    if submitted and email and api_key and query:
        searcher = PubMedSearch(email=email, api_key=api_key)
        results = searcher.search_articles(query, max_results, start_year, end_year)

        if results:
            # Converte para DataFrame
            df = pd.DataFrame(results)

            # Mostra resultados em cards
            for _, article in df.iterrows():
                with st.container():
                    st.markdown(f"""
                    ### {article['title']}
                    **Autores:** {article['authors']}  
                    **Journal:** {article['journal']} ({article['year']})  
                    **Abstract:** {article['abstract']}
                    
                    [Ver no PubMed]({article['url']}) | DOI: {article['doi']}
                    """)
                    st.markdown("---")

            # Opção para download
            st.download_button(
                "📥 Download resultados (CSV)",
                df.to_csv(index=False).encode('utf-8'),
                f"pubmed_results_{int(time.time())}.csv",
                "text/csv",
                key='download-csv'
            )
        else:
            st.warning("Nenhum resultado encontrado.")

    elif submitted:
        st.error("Por favor, preencha email, API token e termo de busca.")

if __name__ == "__main__":
    main()