<template>
  <div>
    <template v-if="displayedPosts.length">
      <PostCard
        v-for="(post, index) in displayedPosts"
        :key="post.id"
        :post="post"
        :index="index"
        :ordinal="(currentPage - 1) * props.perPage + index + 1"
        :category-label="formatCategory(post.category)"
        :article-path="getArticlePath(post)"
        @tag-click="goTag"
      />
    </template>

    <EmptyState v-else icon="fas fa-inbox" show-home>
      {{ t('noPosts') }}
    </EmptyState>

    <Pagination
      :current-page="currentPage"
      :total-pages="totalPages"
      :middle-pages="middlePages"
      :show-first-page="showFirstPage"
      :show-last-page="showLastPage"
      :show-first-ellipsis="showFirstEllipsis"
      :show-last-ellipsis="showLastEllipsis"
      :on-go-to-page="goToPage"
      :on-prev="prevPage"
      :on-next="nextPage"
    />
  </div>
</template>

<script setup lang="ts">
import { computed, onBeforeUnmount, onMounted, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import EmptyState from '@/components/common/EmptyState.vue'
import Pagination from '@/components/common/Pagination.vue'
import PostCard from '@/components/common/PostCard.vue'
import { useLocalizedContent } from '@/composables/useLocalizedContent'
import { useRoutePagination } from '@/composables/useRoutePagination'
import { useTagNavigation } from '@/composables/useTagNavigation'
import { loadNotes, loadCategories } from '@/utils/contentLoader'
import { toArticle } from '@/utils/navigation'
import { getVisiblePages } from '@/utils/pagination'
import type { Post } from '@/types'

const props = withDefaults(defineProps<{ docs: Post[]; perPage?: number }>(), { perPage: 6 })

const { t, locale } = useI18n()
const { goTag } = useTagNavigation()

const maxVisiblePages = ref(5)
const { data: notesFlat } = useLocalizedContent(() => loadNotes(), [])
const { data: categoriesData } = useLocalizedContent(() => loadCategories(), [])

const totalPages = computed(() => Math.max(1, Math.ceil(props.docs.length / props.perPage)))
const { currentPage, goToPage, prevPage, nextPage } = useRoutePagination(
  computed(() => props.docs),
  () => totalPages.value,
)

const displayedPosts = computed(() => {
  const start = (currentPage.value - 1) * props.perPage
  return props.docs.slice(start, start + props.perPage)
})

const allVisiblePages = computed(() =>
  getVisiblePages(currentPage.value, totalPages.value, maxVisiblePages.value),
)

const middlePages = computed(() =>
  allVisiblePages.value.filter(
    (page): page is number => typeof page === 'number' && page !== 1 && page !== totalPages.value,
  ),
)

const showFirstPage = computed(() => allVisiblePages.value.includes(1))
const showLastPage = computed(() => allVisiblePages.value.includes(totalPages.value))
const showFirstEllipsis = computed(
  () => allVisiblePages.value[0] === '...' || allVisiblePages.value[1] === '...',
)
const showLastEllipsis = computed(
  () => allVisiblePages.value[allVisiblePages.value.length - 2] === '...',
)

const categoryTitleMap = computed(() => {
  const map: Record<string, string> = {}
  categoriesData.value.forEach((section) =>
    section.items.forEach((item) =>
      item.categories.forEach((cat) => {
        if (cat.key && cat.title) map[cat.key] = cat.title
      }),
    ),
  )
  return map
})

function getArticlePath(post: Post) {
  const found = notesFlat.value.find((item) => item.title === post.title)
  if (found?.relativePath) return toArticle(`notes/${found.relativePath}.md`)

  const basePath = post.category[1] || 'notes'
  const fileName = post.title.toLowerCase().replace(/[^a-z0-9]/g, '-')
  let articlePath = `${basePath}/${fileName}.md`
  if (locale.value === 'en-US') articlePath = articlePath.replace('.md', '-en.md')
  return toArticle(articlePath)
}

function formatCategory(cat: Post['category']) {
  const [top, subKey] = cat
  if (!subKey) return top
  const subTitle = categoryTitleMap.value[subKey] || subKey
  return `${top} / ${subTitle}`
}

function handleResize() {
  maxVisiblePages.value = window.innerWidth < 480 ? 3 : 5
}

onMounted(() => {
  handleResize()
  window.addEventListener('resize', handleResize)
})

onBeforeUnmount(() => {
  window.removeEventListener('resize', handleResize)
})
</script>
