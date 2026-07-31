<!-- SearchModal.vue -->
<template>
  <Teleport to="body">
    <div class="search-overlay" @click.self="close" @keydown.esc="close">
      <div class="search-panel" role="dialog" aria-modal="true" :aria-label="t('search')">
        <div class="search-input-row d-flex align-items-center">
          <i class="fas fa-search search-icon"></i>
          <input ref="searchInput" v-model="keyword" class="search-input" type="text" inputmode="search"
            :placeholder="t('searchPlaceholder')" @keydown.enter.prevent="openFirstResult" />
          <button class="search-close btn-icon" :aria-label="t('close')" @click="close">
            <span class="search-close-icon" aria-hidden="true"></span>
          </button>
        </div>

        <div v-if="errorMsg" class="search-status">{{ errorMsg }}</div>
        <div v-else-if="keyword && !searching && results.length === 0" class="search-status">
          {{ t('searchNoResults') }}
        </div>

        <ul v-if="results.length" class="search-results list-unstyled mb-0">
          <li v-for="item in results" :key="item.id" class="search-result-item">
            <button class="search-result-btn w-100 text-start" @click="goToResult(item)">
              <div class="d-flex justify-content-between align-items-center gap-2">
                <span class="search-result-title">{{ item.title }}</span>
                <span v-if="item.tags && item.tags.length" class="search-result-path">
                  {{ item.tags.slice(0, 2).join(' / ') }}
                </span>
              </div>
              <p class="search-result-snippet mb-0">{{ snippet(item) }}</p>
            </button>
          </li>
        </ul>
      </div>
    </div>
  </Teleport>
</template>

<script setup lang="ts">
import { nextTick, ref, watch } from 'vue'
import { useI18n } from 'vue-i18n'
import { useRouter } from 'vue-router'
import MiniSearch from 'minisearch'
import { toLocalePath } from '@/utils/localePath'

interface SearchDoc {
  id: string
  title: string
  tags: string[]
  path: string
  description: string
  content: string
}

const emit = defineEmits(['close'])

const { t, locale } = useI18n()
const router = useRouter()

const keyword = ref('')
const results = ref<SearchDoc[]>([])
const searching = ref(false)
const errorMsg = ref('')
const searchInput = ref<HTMLInputElement | null>(null)

let engine: MiniSearch<SearchDoc> | null = null

const cjkTokenize = (text: string) => {
  const words = text.toLowerCase().split(/[^a-z0-9\u4e00-\u9fa5]+/).filter(Boolean)
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
  const mods = import.meta.glob('/src/content/search-index*.json')
  const key = locale.value === 'en-US' ? '/src/content/search-index-en.json' : '/src/content/search-index.json'
  const loader = mods[key]
  if (!loader) return
  const mod = await loader()
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
  const kw = keyword.value.trim()
  if (!kw) {
    results.value = []
    return
  }
  searching.value = true
  errorMsg.value = ''
  try {
    await ensureEngine()
    if (!engine) {
      errorMsg.value = t('searchUnavailable')
      results.value = []
      return
    }
    results.value = engine.search(kw) as unknown as SearchDoc[]
  } catch (e) {
    console.error('Search failed:', e)
    errorMsg.value = t('searchUnavailable')
    results.value = []
  } finally {
    searching.value = false
  }
}

function snippet(item: SearchDoc) {
  const text = item.description || item.content || item.title
  return text.length > 140 ? `${text.slice(0, 140)}...` : text
}

function goToResult(item: SearchDoc) {
  close()
  router.push(toLocalePath(`/article/${item.path}`))
}

function openFirstResult() {
  if (results.value.length > 0) goToResult(results.value[0])
}

function close() {
  emit('close')
}

watch(keyword, () => runSearch())

nextTick(() => {
  searchInput.value?.focus()
})
</script>

<style scoped>
.search-overlay {
  position: fixed;
  inset: 0;
  z-index: 1080;
  background: var(--app-backdrop-bg);
  display: flex;
  align-items: flex-start;
  justify-content: center;
  padding-top: 12vh;
}

.search-panel {
  width: min(640px, 92vw);
  background: var(--app-card-bg);
  border-radius: 0.75rem;
  box-shadow: var(--app-offcanvas-shadow);
  padding: 1rem;
  max-height: 70vh;
  overflow-y: auto;
}

.search-input-row {
  gap: 0.75rem;
}

.search-icon {
  color: var(--app-text-muted);
}

.search-input {
  flex: 1;
  border: none;
  outline: none;
  background: transparent;
  font-size: 1.05rem;
  color: var(--app-text);
}

.search-close {
  border: none;
  background: transparent;
  padding: 0.25rem;
  border-radius: 0.375rem;
  line-height: 1;
  color: var(--app-text-muted);
  transition: color 0.15s ease-in-out;
}

.search-close:hover {
  color: var(--app-text);
}

.search-close-icon {
  position: relative;
  display: inline-block;
  width: 14px;
  height: 14px;
}

.search-close-icon::before,
.search-close-icon::after {
  content: '';
  position: absolute;
  left: 50%;
  top: 50%;
  width: 100%;
  height: 1.5px;
  border-radius: 2px;
  background: currentColor;
}

.search-close-icon::before {
  transform: translate(-50%, -50%) rotate(45deg);
}

.search-close-icon::after {
  transform: translate(-50%, -50%) rotate(-45deg);
}

.search-status {
  padding: 1rem 0.25rem 0.25rem;
  color: var(--app-text-muted);
  font-size: 0.95rem;
}

.search-results {
  margin-top: 0.75rem;
}

.search-result-item + .search-result-item {
  border-top: 1px solid var(--app-border);
}

.search-result-btn {
  background: transparent;
  border: none;
  padding: 0.6rem 0.25rem;
  border-radius: 0.5rem;
  color: var(--app-text);
  transition: background-color 0.15s ease;
}

.search-result-btn:hover,
.search-result-btn:focus-visible {
  background-color: var(--app-primary-bg-subtle);
  outline: none;
}

.search-result-title {
  font-weight: 600;
  font-size: 1rem;
  color: var(--app-primary);
}

.search-result-path {
  font-size: 0.8rem;
  color: var(--app-text-muted);
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
}

.search-result-snippet {
  margin-top: 0.25rem;
  font-size: 0.9rem;
  color: var(--app-text-secondary);
  display: -webkit-box;
  -webkit-line-clamp: 2;
  -webkit-box-orient: vertical;
  overflow: hidden;
}
</style>
