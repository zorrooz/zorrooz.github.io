<!-- PostList.vue -->
<template>
  <div>
    <div v-for="post in displayedPosts" :key="post.id" class="col-12 mb-3">
      <div class="card shadow-sm border-0" :style="{ backgroundColor: 'var(--app-card-bg)' }">
        <div class="card-body p-4 d-flex flex-column flex-md-row gap-4">
          <div class="flex-grow-1">
            <small class="meta-text mb-2 d-block" :style="{ color: 'var(--app-text-muted)' }">
              {{ formatDate(post.date) }}
            </small>

            <h5 class="post-title fw-bold mb-2" :style="{ color: 'var(--app-text)' }">
              <router-link :to="getArticlePath(post)" class="text-decoration-none"
                :style="{ color: 'var(--app-text)' }">
                {{ post.title }}
              </router-link>
            </h5>

            <div class="mb-2">
              <small class="meta-text" :style="{ color: 'var(--app-text-muted)' }">
                {{ formatCategory(post.category) }}
              </small>
            </div>

            <p class="mb-3 desc-text">
              {{ post.preview }}
            </p>


            <div class="d-flex flex-wrap gap-2">
              <span v-for="tag in post.tags" :key="tag"
                class="badge tag-badge fw-normal py-1 px-2 rounded-3 cursor-pointer" @click="goTag(tag)">
                # {{ tag }}
              </span>
            </div>
          </div>
        </div>
      </div>
    </div>

    <div class="row mt-4 pb-4" v-if="totalPages > 1">
      <div class="col-12">
        <nav :aria-label="paginationLabel">
          <ul class="pagination justify-content-between align-items-center mb-0">
            <li class="page-item" :class="{ disabled: currentPage === 1 }">
              <button class="page-link d-flex align-items-center border-0 bg-transparent px-3 py-2"
                :style="{ color: 'var(--app-text-muted)' }" @click="prevPage" :disabled="currentPage === 1"
                :aria-label="t('prevPage')">
                <span class="me-1">&lt;</span>
                <span>{{ t('prevPage') }}</span>
              </button>
            </li>

            <li class="page-item">
              <div class="d-flex gap-2">
                <button v-if="showFirstPage" class="page-link d-inline-flex align-items-center justify-content-center"
                  :class="{ 'current-page': currentPage === 1 }" @click="goToPage(1)">
                  1
                </button>

                <span class="d-flex align-items-center" v-if="showFirstEllipsis">...</span>

                <button v-for="page in middlePages" :key="page"
                  class="page-link d-inline-flex align-items-center justify-content-center"
                  :class="{ 'current-page': currentPage === page }" @click="goToPage(page)">
                  {{ page }}
                </button>

                <span class="d-flex align-items-center" v-if="showLastEllipsis">...</span>

                <button v-if="showLastPage && totalPages > 1"
                  class="page-link d-inline-flex align-items-center justify-content-center"
                  :class="{ 'current-page': currentPage === totalPages }" @click="goToPage(totalPages)">
                  {{ totalPages }}
                </button>
              </div>
            </li>

            <li class="page-item" :class="{ disabled: currentPage === totalPages }">
              <button class="page-link d-flex align-items-center border-0 bg-transparent px-3 py-2"
                :style="{ color: 'var(--app-text-muted)' }" @click="nextPage" :disabled="currentPage === totalPages"
                :aria-label="t('nextPage')">
                <span>{{ t('nextPage') }}</span>
                <span class="ms-1">&gt;</span>
              </button>
            </li>
          </ul>
        </nav>
      </div>
    </div>
  </div>
</template>

<script setup lang="ts">
import { computed, nextTick, onBeforeUnmount, onMounted, ref, watch } from 'vue'
import { useI18n } from 'vue-i18n'
import { useRoute, useRouter } from 'vue-router'
import { loadNotes, loadCategories } from '@/utils/contentLoader'

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

function formatDate(dateString: string) {
  const l = locale.value === 'zh-CN' ? 'zh-CN' : 'en-US'
  const options: Intl.DateTimeFormatOptions = { year: 'numeric', month: 'long', day: 'numeric' }
  return new Date(dateString).toLocaleDateString(l, options)
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

  return `/article/${articlePath.replace(/\.md$/, '')}`
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
  router.push({ path: '/', query: q }).catch(() => { })
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
.post-title {
  font-size: 1.4rem;
  font-weight: 700;
}

.desc-text {
  font-size: 1rem;
  line-height: 1.8;
  color: var(--app-text);
}

.meta-text {
  font-size: 0.95rem;
}

.badge {
  font-size: 0.85rem;
  font-weight: 500;
}

.tag-badge {
  color: var(--app-tag-text) !important;
  background-color: var(--app-tag-bg) !important;
}

.page-link {
  transition: all 0.2s ease;
  color: var(--app-pagination-link-text);
  background-color: var(--app-pagination-link-bg);
  border: none;
  outline: none;
  font-weight: 500;
  min-width: 36px;
  height: 36px;
  border-radius: 8px;
  display: flex;
  align-items: center;
  justify-content: center;
}

.page-link:hover {
  background-color: var(--app-primary-bg-subtle);
  color: var(--app-primary);
}

.current-page {
  background-color: var(--app-primary) !important;
  color: var(--app-pagination-current-text) !important;
  font-weight: 600 !important;
}

.page-link:focus {
  box-shadow: none;
  outline: none;
}
</style>
