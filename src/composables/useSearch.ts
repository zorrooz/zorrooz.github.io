import { ref, watch } from 'vue'
import { useI18n } from 'vue-i18n'
import MiniSearch from 'minisearch'

import type { SearchDoc } from '@/types'

export function useSearch() {
  const { t, locale } = useI18n()

  const keyword = ref('')
  const results = ref<SearchDoc[]>([])
  const searching = ref(false)
  const errorMsg = ref('')

  let engine: MiniSearch<SearchDoc> | null = null
  let searchSeq = 0

  const cjkTokenize = (text: string) => {
    const words = text
      .toLowerCase()
      .split(/[^a-z0-9\u4e00-\u9fa5]+/)
      .filter(Boolean)
    const tokens: string[] = []
    for (const word of words) {
      if (word.length > 2 && /[\u4e00-\u9fa5]/.test(word)) {
        for (let i = 0; i <= word.length - 2; i++) tokens.push(word.slice(i, i + 2))
      } else {
        tokens.push(word)
      }
    }
    return tokens
  }

  async function ensureEngine() {
    if (engine) return
    const currentLocale = locale.value
    const mods = import.meta.glob('@data/content/search-index*.json')
    const key = Object.keys(mods).find((k) =>
      k.includes(currentLocale === 'en-US' ? 'search-index-en.json' : 'search-index.json'),
    )
    if (!key) return
    const loader = mods[key]
    if (!loader) return
    const mod = await loader()
    // 加载期间 locale 已切换：丢弃过期索引，由新 locale 的 runSearch 重建
    if (locale.value !== currentLocale) return
    const docs = ((mod as { default?: SearchDoc[] }).default || []) as SearchDoc[]
    engine = new MiniSearch<SearchDoc>({
      fields: ['title', 'description', 'content'],
      storeFields: ['title', 'tags', 'path', 'description'],
      tokenize: cjkTokenize,
      processTerm: (term) => term,
      searchOptions: {
        prefix: true,
        fuzzy: 0.2,
        boost: { title: 3, description: 2, content: 1 },
      },
    })
    engine.addAll(docs)
  }

  async function runSearch() {
    const seq = ++searchSeq
    const kw = keyword.value.trim()
    if (!kw) {
      results.value = []
      return
    }
    searching.value = true
    errorMsg.value = ''
    try {
      await ensureEngine()
      if (seq !== searchSeq) return
      if (!engine) {
        errorMsg.value = t('searchUnavailable')
        results.value = []
        return
      }
      results.value = engine.search(kw) as unknown as SearchDoc[]
    } catch (e) {
      if (seq !== searchSeq) return
      console.error('Search failed:', e)
      errorMsg.value = t('searchUnavailable')
      results.value = []
    } finally {
      if (seq === searchSeq) searching.value = false
    }
  }

  watch(keyword, () => runSearch())
  watch(locale, () => {
    engine = null
    runSearch()
  })

  return { keyword, results, searching, errorMsg, runSearch }
}
