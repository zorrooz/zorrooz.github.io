<!-- PostList.vue -->
<template>
  <div>
    <article v-for="(post, index) in displayedPosts" :key="post.id" class="post-item" v-reveal
      :style="{ '--reveal-delay': Math.min(index, 5) * 40 + 'ms' }">
      <span class="post-item__index num" aria-hidden="true">{{ String((currentPage - 1) * props.perPage + index + 1).padStart(2, '0') }}</span>

      <router-link :to="getArticlePath(post)" class="post-item__main">
        <div class="post-item__meta">
          <span class="post-item__cat">{{ formatCategory(post.category) }}</span>
          <span class="divider-v"></span>
          <span class="post-item__date"><i class="fas fa-calendar-alt"></i>{{ post.date }}</span>
          <span class="divider-v"></span>
          <span class="post-item__words"><i class="fas fa-file-lines"></i>{{ post.wordCount }} {{ t('wordUnit') }}</span>
          <span class="divider-v"></span>
          <span class="post-item__reading"><i class="fas fa-clock"></i>{{ readingTime(post.wordCount) }}</span>
        </div>
        <h3 class="post-item__title">{{ post.title }}</h3>
        <p class="post-item__preview">{{ post.preview }}</p>
      </router-link>

      <div class="post-item__tags">
        <span v-for="tag in post.tags" :key="tag" class="post-item__tag" @click="goTag(tag)">{{ tag }}</span>
      </div>

      <router-link :to="getArticlePath(post)" class="post-item__arrow" :aria-label="post.title">
        <i class="fas fa-arrow-right"></i>
      </router-link>
    </article>

    <nav class="pagination" v-if="totalPages > 1" :aria-label="paginationLabel">
      <div class="pagination__side">
        <button v-if="currentPage > 1" class="page-btn page-btn--nav" @click="prevPage"
          :aria-label="t('prevPage')">
          <i class="fas fa-chevron-left"></i>
        </button>
      </div>

      <div class="pagination__pages">
        <button v-if="showFirstPage" class="page-btn" :class="{ 'page-btn--active': currentPage === 1 }"
          @click="goToPage(1)">1</button>
        <span class="page-ellipsis" v-if="showFirstEllipsis">...</span>

        <button v-for="page in middlePages" :key="page" class="page-btn"
          :class="{ 'page-btn--active': currentPage === page }" @click="goToPage(page)">
          {{ page }}
        </button>

        <span class="page-ellipsis" v-if="showLastEllipsis">...</span>
        <button v-if="showLastPage && totalPages > 1" class="page-btn"
          :class="{ 'page-btn--active': currentPage === totalPages }" @click="goToPage(totalPages)">
          {{ totalPages }}
        </button>
      </div>

      <div class="pagination__side pagination__side--right">
        <button v-if="currentPage < totalPages" class="page-btn page-btn--nav" @click="nextPage"
          :aria-label="t('nextPage')">
          <i class="fas fa-chevron-right"></i>
        </button>
      </div>
    </nav>
  </div>
</template>

<script setup lang="ts">
import { computed, nextTick, onBeforeUnmount, onMounted, ref, watch } from 'vue'
import { useI18n } from 'vue-i18n'
import { useRoute, useRouter } from 'vue-router'
import { loadNotes, loadCategories } from '@/utils/contentLoader'
import { toLocalePath } from '@/utils/localePath'
import { readingTimeMinutes } from '@/utils/readingTime'

const props = withDefaults(defineProps<{ docs: any[]; perPage?: number }>(), { perPage: 6 })

const { t, locale } = useI18n()
const route = useRoute()
const router = useRouter()

const currentPage = ref(1)
const maxVisiblePages = ref(5)
const notesFlat = ref<any[]>([])
const categoriesData = ref<any[]>([])

const paginationLabel = computed(() => t('pagination'))

const totalPages = computed(() => Math.max(1, Math.ceil(props.docs.length / props.perPage)))

const displayedPosts = computed(() => {
  const start = (currentPage.value - 1) * props.perPage
  const end = start + props.perPage
  return props.docs.slice(start, end)
})

const allVisiblePages = computed(() => {
  const pages: number[] = []
  const total = totalPages.value
  const current = currentPage.value
  const maxShow = maxVisiblePages.value

  if (total <= maxShow) {
    for (let i = 1; i <= total; i++) {
      pages.push(i)
    }
    return pages
  }

  pages.push(current)

  let i = 1
  while (pages.length < maxShow) {
    if (current - i >= 1 && !pages.includes(current - i)) {
      pages.push(current - i)
    }
    if (pages.length >= maxShow) break

    if (current + i <= total && !pages.includes(current + i)) {
      pages.push(current + i)
    }
    i++
  }

  if (!pages.includes(1)) {
    pages.pop()
    pages.push(1)
  }
  if (pages.length < maxShow && !pages.includes(total)) {
    pages.push(total)
  } else if (pages.length >= maxShow && !pages.includes(total)) {
    pages.pop()
    pages.push(total)
  }

  return pages.sort((a, b) => a - b)
})

const middlePages = computed(() =>
  allVisiblePages.value.filter(page => page !== 1 && page !== totalPages.value)
)

const showFirstPage = computed(() => allVisiblePages.value.includes(1))
const showLastPage = computed(() => allVisiblePages.value.includes(totalPages.value))
const showFirstEllipsis = computed(() =>
  totalPages.value > maxVisiblePages.value && allVisiblePages.value[1] > 2
)
const showLastEllipsis = computed(() =>
  totalPages.value > maxVisiblePages.value &&
  allVisiblePages.value[allVisiblePages.value.length - 2] < totalPages.value - 1
)

const categoryTitleMap = computed(() => {
  const map: Record<string, string> = {}
  try {
    (categoriesData.value || []).forEach((section: any) => {
      (section.items || []).forEach((item: any) => {
        (item.categories || []).forEach((cat: any) => {
          if (cat && cat.key && cat.title) map[cat.key] = cat.title
        })
      })
    })
  } catch (e) {
    console.error(e)
  }
  return map
})

function loadData() {
  try {
    notesFlat.value = loadNotes() || []
    categoriesData.value = loadCategories() || []
  } catch (error) {
    console.error('Failed to load data:', error)
    notesFlat.value = []
    categoriesData.value = []
  }
}

function getArticlePath(post: any) {
  let articlePath = ''

  const found = Array.isArray(notesFlat.value) ? notesFlat.value.find(item => item.title === post.title) : null
  if (found && found.relativePath) {
    articlePath = `notes/${found.relativePath}.md`
  } else {
    const isEnglish = locale.value === 'en-US'

    const basePath = post.category?.[1] || 'notes'
    const fileName = post.title.toLowerCase().replace(/[^a-z0-9]/g, '-')
    articlePath = `${basePath}/${fileName}.md`

    if (isEnglish) {
      articlePath = articlePath.replace('.md', '-en.md')
    }
  }

  return toLocalePath(`/article/${articlePath.replace(/\.md$/, '')}`)
}

function goToPage(page: number) {
  if (page >= 1 && page <= totalPages.value) {
    currentPage.value = page
    const q = { ...route.query, page: String(page) }
    router.push({ path: route.path, query: q }).catch(() => { })
    nextTick(() => window.scrollTo({ top: 0, behavior: 'smooth' }))
  }
}

function prevPage() {
  if (currentPage.value > 1) {
    const page = currentPage.value - 1
    currentPage.value = page
    const q = { ...route.query, page: String(page) }
    router.push({ path: route.path, query: q }).catch(() => { })
    nextTick(() => window.scrollTo({ top: 0, behavior: 'smooth' }))
  }
}

function nextPage() {
  if (currentPage.value < totalPages.value) {
    const page = currentPage.value + 1
    currentPage.value = page
    const q = { ...route.query, page: String(page) }
    router.push({ path: route.path, query: q }).catch(() => { })
    nextTick(() => window.scrollTo({ top: 0, behavior: 'smooth' }))
  }
}

function goTag(tag: string) {
  if (!tag) return
  const q = { ...route.query, tag: tag, page: '1' }
  const prefix = locale.value === 'en-US' ? '/en' : '/zh'
  router.push({ path: `${prefix}/`, query: q }).catch(() => { })
  nextTick(() => window.scrollTo({ top: 0, behavior: 'smooth' }))
}

function handleResize() {
  maxVisiblePages.value = window.innerWidth < 480 ? 3 : 5
}

function formatCategory(catArr: any) {
  if (!Array.isArray(catArr) || catArr.length === 0) return ''
  const top = catArr[0] || ''
  const subKey = catArr[1]
  if (!subKey) return top
  const subTitle = categoryTitleMap.value[subKey] || subKey
  return `${top} / ${subTitle}`
}

function readingTime(wordCount: number) {
  return t('postReadingTime', { minutes: readingTimeMinutes(wordCount) })
}

loadData()

watch(() => props.docs, () => { currentPage.value = 1 })

watch(() => route.query.page, (newVal) => {
  const p = parseInt(String(newVal))
  const page = Number.isFinite(p) && p >= 1 ? Math.min(p, totalPages.value) : 1
  if (page !== currentPage.value) {
    currentPage.value = page
    nextTick(() => window.scrollTo({ top: 0, behavior: 'smooth' }))
  }
})

watch(locale, () => loadData())

onMounted(() => {
  const p = parseInt(String(route.query.page))
  currentPage.value = Number.isFinite(p) && p >= 1 ? Math.min(p, totalPages.value) : 1
  handleResize()
  window.addEventListener('resize', handleResize)
})

onBeforeUnmount(() => {
  window.removeEventListener('resize', handleResize)
})
</script>

<style scoped>
.post-item {
  display: grid;
  grid-template-columns: 64px 1fr auto;
  grid-template-areas:
    'index main arrow'
    'index tags arrow';
  align-items: center;
  gap: var(--sp-2) var(--sp-4);
  padding: var(--sp-6) 0;
  border-bottom: 1px solid var(--line);
  position: relative;
}

.post-item:last-child {
  border-bottom: none;
}

.post-item:hover {
  background: transparent;
}

.post-item__index {
  grid-area: index;
  align-self: start;
  padding-top: 6px;
  font-size: 22px;
  font-weight: 700;
  color: var(--line-strong);
  letter-spacing: 0.02em;
  font-variant-numeric: tabular-nums;
  transition: color 0.14s ease;
}

.post-item__main {
  grid-area: main;
  min-width: 0;
  display: block;
  text-decoration: none;
  color: inherit;
}

.post-item__meta {
  display: flex;
  align-items: center;
  gap: var(--sp-3);
  flex-wrap: wrap;
  font-size: var(--text-xs);
  color: var(--fg-3);
}

.post-item__meta i {
  font-size: 12px;
  color: var(--primary-muted);
  margin-right: 6px;
  transition: color 0.14s ease;
}

.post-item:hover .post-item__meta i {
  color: var(--primary);
}

.post-item__cat {
  font-weight: 600;
  color: var(--fg-2);
}

.post-item__date {
  letter-spacing: 0.01em;
  font-size: var(--text-xs);
}

.post-item__words {
  color: var(--fg-3);
}

.post-item__title {
  font-size: var(--text-xl);
  font-weight: 600;
  letter-spacing: -0.015em;
  line-height: 1.35;
  margin: var(--sp-3) 0 var(--sp-2);
  color: var(--fg);
  transition: color 0.14s ease;
}

.post-item:hover .post-item__title {
  color: var(--primary);
}

.post-item__preview {
  font-size: var(--text-md);
  color: var(--fg-2);
  line-height: 1.65;
  margin: 0;
  display: -webkit-box;
  -webkit-line-clamp: 2;
  -webkit-box-orient: vertical;
  overflow: hidden;
  max-width: 640px;
}

.post-item__tags {
  grid-area: tags;
  display: flex;
  gap: var(--sp-2);
  margin-top: var(--sp-3);
  flex-wrap: wrap;
}

.post-item__tag {
  font-size: var(--text-xs);
  font-weight: 500;
  padding: 3px 11px;
  border-radius: var(--radius-pill);
  background: var(--surface-2);
  color: var(--fg-2);
  cursor: pointer;
  transition: color 0.14s ease, background-color 0.14s ease;
}

.post-item__tag:hover {
  color: var(--primary);
  background: var(--tint);
}

.post-item__arrow {
  grid-area: arrow;
  align-self: center;
  display: flex;
  align-items: center;
  justify-content: center;
  color: var(--fg-3);
  font-size: 15px;
  text-decoration: none;
  padding: 6px;
  opacity: 0;
  transition: opacity var(--dur-base) ease, color var(--dur-fast) ease;
}

.post-item:hover .post-item__arrow {
  opacity: 1;
}

.post-item__arrow:hover {
  color: var(--primary);
}

.pagination {
  display: grid;
  grid-template-columns: 1fr auto 1fr;
  align-items: center;
  gap: var(--sp-2);
  padding: var(--sp-12) 0 var(--sp-4);
}

.pagination__side {
  justify-self: start;
}

.pagination__side--right {
  justify-self: end;
}

.pagination__pages {
  display: flex;
  align-items: center;
  gap: var(--sp-2);
  justify-self: center;
}

.page-btn--nav {
  min-width: 34px;
  height: 34px;
  border: none;
  background: transparent;
  color: var(--fg-3);
  border-radius: 50%;
}

.page-btn--nav:hover {
  border-color: transparent;
  color: var(--primary);
}

.page-btn:disabled {
  opacity: 0.4;
  cursor: not-allowed;
}

.page-ellipsis {
  font-size: 14px;
  color: var(--fg-3);
  padding: 0 4px;
}

@media (max-width: 767px) {
  .post-item {
    padding: var(--sp-5) 0;
    grid-template-columns: 1fr auto;
    grid-template-areas:
      'main arrow'
      'tags arrow';
  }

  .post-item__index {
    display: none;
  }

  .post-item__title {
    font-size: 17px;
  }

  .post-item__arrow {
    opacity: 1;
  }
}
</style>
